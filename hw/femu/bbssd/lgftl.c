#include "lgftl.h"

#include <math.h>

static void *ftl_thread(void *arg);

/* --- linear regression functions ---*/

typedef struct linear_model_t {
    float w;
    float b;
} linear_model_t;

typedef struct ransac_t {
    uint32_t samplesize;
    float maxerror;
    
    uint32_t bestfit_cnt;
    linear_model_t bestmodel;
} ransac_t;

static void LinearRegression(uint32_t *x, uint32_t *y, int num, linear_model_t *lm) {  
        float xsum=0, ysum=0, xxsum=0, xysum=0;  
        for(int i = 0; i < num; i++) {  
            xxsum += x[i]*x[i];  
            xsum += x[i];  
            xysum += x[i]*y[i];  
            ysum += y[i];  
        }  
        lm->w = (xysum * num - xsum * ysum) / (xxsum * num - xsum * xsum);   
        lm->b = (ysum - lm->w * xsum) / num;  
}  

static void EvaluateError(uint32_t *x, uint32_t *y, int num, linear_model_t *lm, float *error) {
    for (int i = 0; i < num; i++) {
        float estimated_y = lm->w * x[i] + lm->b;
        error[i] = estimated_y - y[i];
    }
}

static uint32_t CntInliers(float *error, uint32_t num, float maxerror) {
    uint32_t inliers = 0;
    for (uint32_t i = 0; i < num; i++) {
        inliers += (fabsf(error[i]) < maxerror ? 1 : 0);
    }
    return inliers;
}

static uint32_t IterationNum(float iratio, uint32_t samplesize, float successrate) {
    return log(1 - successrate) / log(1 - pow(iratio, samplesize));
}

static void SampleTrainingSet(uint32_t *x, uint32_t *y, uint32_t num, uint32_t samplesize, uint32_t *s_x, uint32_t *s_y) {
    srand(time(NULL));
    for (uint32_t s = 0; s < samplesize; s++) {
        uint32_t r = rand() % num;
        s_x[s] = x[r];
        s_y[s] = y[r];
    }
}

static void RANSAC_fit (uint32_t *x, uint32_t *y, int num, ransac_t *model) {
    /* init */
    model->bestfit_cnt = 0;
    uint32_t *sample_x = g_malloc(sizeof(uint32_t) * model->samplesize); 
    uint32_t *sample_y = g_malloc(sizeof(uint32_t) * model->samplesize);
    float *error = g_malloc(sizeof(float) * model->samplesize);
    
    uint32_t numiteration = IterationNum(0.7, model->samplesize, 0.99);
    uint32_t iter = 0;
    while (iter++ < numiteration) {
        // create smapled data
        SampleTrainingSet(x, y, num, model->samplesize, sample_x, sample_y);
        // fit the modle on sampled data
        linear_model_t tmp_lm;
        LinearRegression(sample_x, sample_y, model->samplesize, &tmp_lm);
        // evaluate the model
        EvaluateError(x, y, num, &tmp_lm, error);
        uint32_t inliers = CntInliers(error, num, model->maxerror);
        if (inliers > model->bestfit_cnt) {
            model->bestmodel = tmp_lm;
        }
    }

    // final round to fit a model for all inliers
    uint32_t *inliers_x = g_malloc(sizeof(uint32_t) * model->bestfit_cnt);
    uint32_t *inliers_y = g_malloc(sizeof(uint32_t) * model->bestfit_cnt);
    uint32_t i = 0;
    for (uint32_t j = 0; j < num; j++) {
        if (error[j] < model->maxerror) {
            inliers_x[i] = error[j];
            inliers_y[i] = error[j];
            if (i++ == model->bestfit_cnt) {
                break;
            }
        }
    }
    LinearRegression(inliers_x, inliers_y, i, &model->bestmodel);
    
    g_free(sample_x);
    g_free(sample_y);
    g_free(error);
    g_free(inliers_x);
    g_free(inliers_y);
}

/* --- Quick Sort functions */

void Exchange(uint32_t *a, uint32_t *b){
    uint32_t tmp = *a;
    *a = *b;
    *b = tmp;
}

void QsortKV(uint32_t *keys, uint32_t *value, int left, int right){  
    if(left >= right) return; //終止條件
    int l = left + 1; //左
    int r = right; //右
    int key = keys[left];
    while(1) {
        while( l <= right){
            if(keys[l] > key) break;
            l++;
        }
        while(r>left){
            if(keys[r] < key) break;
            r--;
        }
        if(l>r) break;
        Exchange(&keys[l], &keys[r]);
        Exchange(&value[l], &value[r]);
    }
    //key 和 相遇的值 互換
    Exchange(&keys[left], &keys[r]);
    Exchange(&value[left], &value[r]);
    //分小組繼續進行
    QSORT(keys, value, left, r-1);
    QSORT(keys, value, r+1, right);
}

/* --- process hash --- */
static inline uint32_t cmt_hash(uint64_t lpn) {
    return lpn % CMT_HASH_SIZE;
}

static inline uint32_t tp_hash(uint64_t tvpn) {
    return tvpn % TP_HASH_SIZE;
}

static inline uint32_t outlier_hash(uint16_t lpn_offset) {
    return lpn_offset % OUTLIER_HASH_SIZE;
}

/* --- mapping table related functions --- */
static inline uint32_t subspace_idx (uint64_t lpn);
static inline uint64_t subspace_baselpn (uint64_t lpn);
static uint64_t virtual_tpage_num(struct gtd_entry *gtd, uint64_t lpn, uint32_t subspace_boundary_lower);
static uint32_t get_covered_tpage_len(struct gtd_entry *gtd, uint64_t tvpn, uint32_t subspace_boundary_upper);
static struct cmt_entry *find_cmt_entry(struct cmt_mgmt *cm, uint64_t lpn);
static struct TPnode *find_tpnode_by_tvpn(struct cmt_mgmt *cm, uint64_t tvpn);
static struct TPnode *find_tpnode_by_lpn(struct cmt_mgmt *cm, struct gtd_entry *gtd, uint64_t lpn);
static void cmt_evict(struct cmt_mgmt *cm);
static cmt_entry *allocate_cmt_entry(struct cmt_mgmt *cm);
static void insert_tp_hashtable(struct cmt_mgmt *cm, TPnode *tpnode);
static void reclaim_tpnode(struct cmt_mgmt *cm, TPnode *tpnode);
static struct TPnode  *allocate_TPnode(struct cmt_mgmt *cm);
static void insert_cmt_entry(struct cmt_mgmt *cm, struct gtd_entry *gtd, cmt_entry *mapping_entry);
static void insert_lm_entry(struct cmt_mgmt *cm, struct gtd_entry *gtd, cmt_entry *mapping_entry);


static inline uint32_t subspace_idx (uint64_t lpn) {
    uint32_t sub_space_idx = lpn / SUB_SPACE_SIZE;
    return sub_space_idx;
}

static inline uint64_t subspace_baselpn (uint64_t lpn) {
    uint64_t base_lpn = (lpn) * SUB_SPACE_SIZE;
    return base_lpn;
}

/* subspace_boundary_lower is the inclusive lower boundary */
static uint64_t virtual_tpage_num(struct gtd_entry *gtd, uint64_t lpn, uint32_t subspace_boundary_lower) {
    uint64_t gtd_offset = lpn / ENT_PER_TP;
    while (gtd_offset > subspace_boundary_lower && gtd[gtd_offset].translation_ppa.vppa == gtd[gtd_offset-1].translation_ppa.vppa) {
        gtd_offset--;
    }
    return gtd_offset;
}

/* subspace_boundary_upper is the exclusive upper boundary */
static uint32_t get_covered_tpage_len(struct gtd_entry *gtd, uint64_t tvpn, uint32_t subspace_boundary_upper) {
    uint32_t covered_gtd_len = 1;
    while (gtd[tvpn].translation_ppa.vppa == gtd[tvpn++].translation_ppa.vppa && tvpn < subspace_boundary_upper) {
        covered_gtd_len++;
    }
    return covered_gtd_len;
}

static struct cmt_entry *find_cmt_entry(struct cmt_mgmt *cm, uint64_t lpn) {
    uint32_t pos = cmt_hash(lpn);
    cmt_entry *mapping_entry = NULL;
    QTAILQ_FOREACH(mapping_entry, &cm->hash_mapping_table[pos], h_entry) {
        bool is_found = false;
        switch(mapping_entry->type){
            case SINGLE_MAPPING:
                if (mapping_entry->single.lpn == lpn) {
                    is_found = true;
                }
                break;
            case LINEAR_MODEL:
                break;
        }
        if (is_found) {
            break;
        }
    }
    return mapping_entry;
}

static struct TPnode *find_tpnode_by_tvpn(struct cmt_mgmt *cm, uint64_t tvpn) {
    uint32_t pos = tp_hash(tvpn);
    TPnode *tpnode = NULL;
    QTAILQ_FOREACH(tpnode, &cm->hash_tp_table[pos], h_entry) {
        if (tpnode->tvpn == tvpn) {
            break;
        }
    }
    return tpnode;
}

static struct TPnode *find_tpnode_by_lpn(struct cmt_mgmt *cm, struct gtd_entry *gtd, uint64_t lpn) {
    uint64_t subspace_base_lpn = subspace_baselpn(lpn);
    uint32_t subspace_boundary_lower = subspace_base_lpn / ENT_PER_TP;
    uint64_t tvpn = virtual_tpage_num(gtd, lpn, subspace_boundary_lower);
    return find_tpnode_by_tvpn(cm, tvpn);
}

static struct cmt_entry *find_lm_entry(struct cmt_mgmt *cm, struct gtd_entry *gtd, uint64_t lpn) {
    uint32_t offset_in_subspace = lpn % SUB_SPACE_SIZE;
    struct TPnode *tpnode = find_tpnode_by_lpn(cm, gtd, lpn);
    struct cmt_entry *lm = NULL;
    if (tpnode) {
        QTAILQ_FOREACH_REVERSE(lm, &tpnode->lm_list, entry) {
            if (lm->lm.start_lpn_offset <= offset_in_subspace) {
                break;
            }
        }
    }
    return lm;
}

static void TPnode_split(struct TPnode *tpnode, struct cmt_mgmt *cm, struct gtd_entry *gtd) {
    if (tpnode->covered_gtd_len == 1) {
        return;
    }
    uint16_t mapping_entry_cnt = 0;
    uint32_t curr_tvpn = tpnode->tvpn, neighbor_tvpn;
    for (int i = 0; i < tpnode->covered_gtd_len; i++) {
        mapping_entry_cnt += (gtd[curr_tvpn+i].lm_cnt + gtd[curr_tvpn+i].outliers_cnt);
    }
    if (mapping_entry_cnt > ENT_PER_MP) {
        // create the right node
        // move the tpnode list
        // move the gtd outliers
        struct TPnode *new_tpnode = allocate_TPnode(cm);
        new_tpnode->tvpn = curr_tvpn + tpnode->covered_gtd_len / 2;
        new_tpnode->covered_gtd_len = tpnode->covered_gtd_len / 2;
        uint64_t lpn_boundary = new_tpnode->tvpn * ENT_PER_TP;
        uint64_t subspace_base_lpn = subspace_baselpn(lpn_boundary);
        struct cmt_entry *curr_ce = NULL, *next_ce = NULL;
        QTAILQ_FOREACH_SAFE(curr_ce, &tpnode->cold_list, entry, next_ce) { 
            if (curr_ce->single.lpn >= lpn_boundary) {
                QTAILQ_REMOVE(&tpnode->cold_list, curr_ce, entry);
                QTAILQ_INSERT_TAIL(&new_tpnode->cold_list, curr_ce, entry);
            }
        }
        bool cutline = false;
        QTAILQ_FOREACH_REVERSE_SAFE(curr_ce, &tpnode->lm_list, entry, next_ce) {
            if (curr_ce->lm.start_lpn_offset + subspace_base_lpn >= lpn_boundary) {
                QTAILQ_REMOVE(&tpnode->lm_list, curr_ce, entry);
                QTAILQ_INSERT_TAIL(&new_tpnode->lm_list, curr_ce, entry);
            } else if (!cutline) { 
                // The first linear model whose start_lpn is less than boundary, its range may cross two TPnode
                // We duplicate it
                cutline = true;
                struct cmt_entry *ce = allocate_cmt_entry(cm);
                
            }
        }
    }
}

static void cmt_evict(struct cmt_mgmt *cm, struct gtd_entry *gtd) {
    /*  find the lru tpnode
        find the clean entry from the tail of the cold list
        if no clean entry, write back the tpage
            if the tpnode is too large, then split
            if the neighbor tpnode can be compacted together, then compact
    */
    struct TPnode *tpnode = QTAILQ_LAST(&cm->TPnode_list);
    struct cmt_entry *ce = NULL;
    bool clean_found = false;
    QTAILQ_FOREACH_REVERSE(ce, &tpnode->cold_list, entry) {
        if (ce->dirty == CLEAN) {
            clean_found = true;
            break;
        }
    }
    if (clean_found) {
        reclaim_cmt_entry(cm, tpnode, ce);
    } else {
        uint16_t mapping_entry_cnt = 0;
        uint32_t curr_tvpn = tpnode->tvpn, neighbor_tvpn;
        for (int i = 0; i < tpnode->covered_gtd_len; i++) {
            mapping_entry_cnt += (gtd[curr_tvpn+i].lm_cnt + gtd[curr_tvpn+i].outliers_cnt);
        }
        if (mapping_entry_cnt > ENT_PER_TP) { // split
            
        }
        if (curr_tvpn % (2 * tpnode->covered_gtd_len)) {
            neighbor_tvpn = curr_tvpn + tpnode->covered_gtd_len;
        }

        
    }
}

static cmt_entry *allocate_cmt_entry(struct cmt_mgmt *cm) {
    // !!! TODO: evict
    struct cmt_entry *ce = g_malloc0(sizeof(struct cmt_entry));
    cm->live_cmt_entry_cnt++;
    cm->free_cmt_entry_cnt--;
    return ce;
}

static void insert_tp_hashtable(struct cmt_mgmt *cm, TPnode *tpnode) {
    uint32_t pos = tp_hash(tpnode->tvpn);
    QTAILQ_INSERT_HEAD(&cm->hash_tp_table[pos], tpnode, h_entry);
}

static void reclaim_tpnode(struct cmt_mgmt *cm, TPnode *tpnode) {
    ftl_assert(tpnode->cmt_entry_cnt == 0 && lm_cnt == 0);
    QTAILQ_REMOVE(&cm->TPnode_list, tpnode, lru_entry); /* detach from the lru list*/
    uint32_t pos = tp_hash(tpnode->tvpn);
    QTAILQ_REMOVE(&cm->hash_tp_table[pos], tpnode, h_entry); /* detach from the hash table */
    g_free(tpnode);
    cm->live_tpnode_cnt--;
    cm->free_cmt_entry_cnt+=2;
}

static struct TPnode  *allocate_TPnode(struct cmt_mgmt *cm) {
    // !!! TODO: evict
    struct TPnode *tpnode = g_malloc0(sizeof(struct TPnode));
    QTAILQ_INIT(&tpnode->lm_list);
    QTAILQ_INIT(&tpnode->cold_list);
    QTAILQ_INIT(&tpnode->hot_list);
    cm->live_tpnode_cnt++;
    cm->free_cmt_entry_cnt-=2;
    return tpnode;
}

static void insert_cmt_entry(struct cmt_mgmt *cm, struct gtd_entry *gtd, cmt_entry *mapping_entry) {
    ftl_assert(mapping_entry->type == SINGLE_MAPPING);
    struct TPnode *tpnode = NULL;
    // insert to the hash table
    uint32_t pos = cmt_hash(mapping_entry->single.lpn);
    QTAILQ_INSERT_HEAD(&cm->hash_mapping_table[pos], mapping_entry, h_entry);
    tpnode = find_tpnode_by_lpn(cm, gtd, mapping_entry->single.lpn);
    // create the tpnode, if not exist
    if (tpnode == NULL) {
        tpnode = allocate_TPnode(cm);
        uint64_t subspace_base_lpn = subspace_baselpn(mapping_entry->single.lpn);
        uint32_t subspace_boundary_lower = subspace_base_lpn / ENT_PER_TP;
        tpnode->tvpn = virtual_tpage_num(gtd, mapping_entry->single.lpn, subspace_boundary_lower);
        uint32_t subspace_boundary_upper = (subspace_base_lpn + SUB_SPACE_SIZE) / ENT_PER_TP;
        tpnode->covered_gtd_len = get_covered_tpage_len(gtd, tpnode->tvpn, subspace_boundary_upper);
    }
    // insert to the TPnode list
    QTAILQ_INSERT_HEAD(&tpnode->cold_list, mapping_entry, entry);
}

static void insert_lm_entry(struct cmt_mgmt *cm, struct gtd_entry *gtd, cmt_entry *mapping_entry) {

}

/*
* Only reclaim cmt_entry from the cold list.
* Before calling this function, the target to be removed should be determined and cleaned.
*/
static void reclaim_cmt_entry(struct cmt_mgmt *cm, TPnode *tpnode, cmt_entry *mapping_entry) {
    QTAILQ_REMOVE(&tpnode->cold_list, mapping_entry, entry); /* detach from the lru list*/
    uint64_t pos = cmt_hash(mapping_entry->single.lpn);
    QTAILQ_REMOVE(&cm->hash_tp_table[pos], mapping_entry, h_entry); /* detach from the hash table */
    cm->live_cmt_entry_cnt--;
    cm->free_cmt_entry_cnt++;
    g_free(mapping_entry);
    if (QTAILQ_EMPTY(&tpnode->lm_list) && QTAILQ_EMPTY(&tpnode->cold_list) && QTAILQ_EMPTY(&tpnode->hot_list)) {
        remove_tpnode(cm, tpnode);
    }
}


/* --- gtd related functions --- */

static void insert_outlier(struct gtd_entry *ge, struct cmt_mgmt *cm, uint16_t lpn_offset) {
    struct outlier *ol = g_malloc(sizeof(struct outlier));
    ol->lpn_offset = lpn_offset;
    uint32_t idx = outlier_hash(lpn_offset);
    QTAILQ_INSERT_HEAD(&ge->outliers_hash_table[idx], ol, h_entry);
    // update the cmt management info
    cm->live_outliers_cnt++;
}

static void remove_outlier(struct gtd_entry *ge, struct cmt_mgmt *cm, uint16_t lpn_offset) {
    uint32_t idx = outlier_hash(lpn_offset);
    struct outlier *ol = NULL, *ol_next = NULL;
    QTAILQ_FOREACH_SAFE(ol, &ge->outliers_hash_table[idx], h_entry, ol_next) {
        if (ol->lpn_offset == lpn_offset) {
            QTAILQ_REMOVE(&ge->outliers_hash_table[idx], ol, h_entry);
            break;
        }
    }
    // update the cmt management info
    cm->live_outliers_cnt--;
}

static bool is_outlier(struct gtd_entry *ge, uint16_t lpn_offset) {
    uint32_t idx = outlier_hash(lpn_offset);
    struct outlier *ol = NULL;
    QTAILQ_FOREACH(ol, &ge->outliers_hash_table[idx], h_entry) {
        if (ol->lpn_offset == lpn_offset) {
            return true;
        }
    }
    return false;
}

/* --- linear group related functions ---*/

struct ppa get_blk(struct ssd *ssd) {
    struct lg_allocate_pointer *lg_ap = &ssd->lgm.lg_ap;
    struct ssd_channel *channels = ssd->ch;
    int ch_num = lg_ap->ch;
    int lun_num = lg_ap->lun;
    int pl_num = lg_ap->pl;
    struct nand_plane *target_plane = &channels[ch_num].lun[lun_num].pl[pl_num];
    struct nand_block *candidate_block = NULL;

    candidate_block = QTAILQ_FIRST(&target_plane->free_blk_list);
    if (!candidate_block) {
        struct ppa noppa = {.ppa = INVALID_PPA};
        return noppa;
    }

    QTAILQ_REMOVE(&target_plane->free_blk_list, candidate_block, entry);  
    candidate_block->is_allocated = true;   
    target_plane->nfreeblks--;
    return candidate_block->blk_ppa;
}

static void advance_lg_allocate_pointer(struct ssd *ssd) {
    struct lg_allocate_pointer *lg_ap = &ssd->lgm.lg_ap;
    struct ssdparams *sp = &ssd->sp;
    lg_ap->ch++;
    if(lg_ap->ch == sp->nchs) {
        lg_ap->ch = 0;
        lg_ap->lun++;
        if (lg_ap->lun == sp->luns_per_ch) {
            lg_ap->lun = 0;
            lg_ap->pl++;
            if (lg_ap->pl == sp->pls_per_lun) {
                lg_ap->pl = 0;
            }
        }
    }
}

struct ppa* allocate_blks(struct ssd *ssd) {
    struct ssdparams *sp = &ssd->sp;
    struct ssd_channel *channels = ssd->ch;
    struct lg_mgmt *lgm= &ssd->lgm;

    int lg_len = sp->blks_per_lg;
    struct ppa *blks = g_malloc0(sizeof(struct ppa) * lg_len);
    for (int i = 0; i < lg_len; i++) {
        // add one blk from current lun to new linear group
        blks[i] = get_blk(ssd);
        while (blks[i].ppa == INVALID_PPA) {
            advance_lg_allocate_pointer(ssd);
            blks[i] = get_blk(ssd);
        }
        advance_lg_allocate_pointer(ssd);
    }
}

struct linear_group* create_lg(struct ssd *ssd, int sub_space_id, int type) {
    struct lg_mgmt *lgm= &ssd->lgm;
    struct ssdparams *sp = &ssd->sp;

    struct linear_group *new_lg = g_malloc0(sizeof(struct linear_group));
    new_lg->is_open = true;
    new_lg->id = lgm->tt_lg;
    new_lg->len = sp->blks_per_lg;
    new_lg->type = type;
    new_lg->blks = allocate_blks(ssd);
    QTAILQ_INIT(&new_lg->reverse_tvpns);

    QTAILQ_INSERT_HEAD(&lgm->lg_list[sub_space_id], new_lg, h_entry);
    lgm->tt_lg++;
    lgm->open_lg_cnt++;
}

void close_lg(struct ssd *ssd, struct linear_group *lg) {
    uint32_t start_offset = lg->start_offset; // inclusive
    uint32_t final_offset = lg->lg_wp.pg * lg->len + lg->lg_wp.blk; // exclusive
    uint64_t base_lpn = lg->subspace_id * SUB_SPACE_SIZE;
    uint32_t num = final_offset - start_offset;
    uint32_t *lpn_offset = g_malloc(sizeof(uint32_t) * num);
    uint32_t *ppn_offset = g_malloc(sizeof(uint32_t) * num);
    bool *precise_bitmap = g_malloc(sizeof(bool) * num);

    // collect training data from the cached mapping table
    struct reverse_tvpn_list *r_tvpn = NULL;
    uint32_t traindata_idx = 0;
    QTAILQ_FOREACH(r_tvpn, &lg->reverse_tvpns, next) {
        struct TPnode *tpnode = find_tpnode_by_tvpn(&ssd->cm, r_tvpn->tvpn);
        struct cmt_entry *mapping_entry = NULL;
        QTAILQ_FOREACH(mapping_entry, &tpnode->hot_list, entry) {
            lpn_offset[traindata_idx] = mapping_entry->single.lpn - base_lpn;
            ppn_offset[traindata_idx] = mapping_entry->single.virtual_ppa.g.lg_offset;
            traindata_idx++;
        }
    }
    ftl_assert(traindata_idx == num - 1);
    QsortKV(lpn_offset, ppn_offset, 0, num);

    // fit the model
    struct ransac_t model = {.bestfit_cnt = 0, .maxerror = 1.0f, .samplesize = num/10};
    RANSAC_fit(lpn_offset, ppn_offset, num, &model);

    linear_model_t *lm = &model.bestmodel;
    uint32_t lm_start_offset = 0, lm_len = 0; 

    /* update the outliers record in GTD
    * two types outliers: (1) outside the model (2) inside the model but unprecise
    */
    uint32_t base_gtd_idx = base_lpn / ENT_PER_TP;
    for (uint32_t i = 0; i < num; i++) {
        uint32_t gtd_idx = lpn_offset[i] / ENT_PER_TP + base_gtd_idx; // calculate the gtd idx for current entry;
        uint16_t gtd_lpn_offset = (lpn_offset[i] + base_lpn) % ENT_PER_TP; // TODO: if base_lpn % ENT_PER_TP == 0, the addition should be reduced
        if (roundf(lm->w * lpn_offset[i] + lm->b) == ppn_offset[i]) { // precise point
            // remove outlier record, if exist
            remove_outlier(&ssd->gtd[gtd_idx], &ssd->cm, gtd_lpn_offset);
        } else {
            // insert outlier record
            insert_outlier(&ssd->gtd[gtd_idx], &ssd->cm, gtd_lpn_offset);
        }
        if (i == 0) {
            // init model para
            lm_start_offset = 0;
        }    
        if (i == num-1 || lpn_offset[i] != lpn_offset[i+1]) {
            lm_len = i - lm_start_offset;
            // segment the model, insert the model into the lm_list
            struct cmt_entry *ce = get_free_cmt_entry(&ssd->cm);
            ce->type = LINEAR_MODEL;
            ce->prefetch = READY;
            ce->lm.slope = lm->w;
            ce->lm.intercept = lm->b;
            ce->lm.lg_id = lg->id;
            ce->lm.start_lpn_offset = lm_start_offset;
            ce->dirty = DIRTY;

            lm_start_offset = i + 1;
        }
        if (lm_len > 0) {

        }
    }
}