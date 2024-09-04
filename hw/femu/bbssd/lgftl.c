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

/* process hash */
static inline uint32_t cmt_hash(uint64_t lpn) {
    return lpn % CMT_HASH_SIZE;
}

static inline uint32_t tp_hash(uint64_t tvpn) {
    return tvpn % TP_HASH_SIZE;
}

static inline uint32_t subspace_idx (uint64_t lpn) {
    uint64_t tvpn = lpn / ENT_PER_TP;
    uint32_t sub_space_idx = tvpn / SUB_SPACE_SIZE;
    return sub_space_idx;
}

static uint64_t virtual_tpage_num (struct gtd_entry *gtd, uint64_t lpn) {
    uint64_t gtd_offset = lpn / ENT_PER_TP;
    while (gtd_offset && gtd[gtd_offset].translation_ppa.vppa == gtd[gtd_offset - 1].translation_ppa.vppa) {
        gtd_offset--;
    }
    return gtd_offset;
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
    uint64_t tvpn = virtual_tpage_num(gtd, lpn);
    return find_tpnode_by_tvpn(cm, tvpn);
}

static struct cmt_entry *find_lm_entry(struct cmt_mgmt *cm, struct gtd_entry *gtd, uint64_t lpn) {
    uint32_t offset_in_subspace = lpn % SUB_SPACE_SIZE;
    struct TPnode *tpnode = find_tpnode_by_lpn(cm, gtd, lpn);
    struct cmt_entry *lm = NULL;
    if (tpnode) {
        QTAILQ_FOREACH(lm, &tpnode->lm_list, entry) {
            if (lm->lm.start_lpn_offset == offset_in_subspace) {
                break;
            }
        }
    }
    return lm;
}

static void insert_cmt_entry(struct cmt_mgmt *cm, cmt_entry *mapping_entry) {
    uint32_t pos = cmt_hash(mapping_entry->single.lpn);
    QTAILQ_INSERT_HEAD(&cm->hash_mapping_table[pos], mapping_entry, h_entry);
}

static void insert_tp_hashtable(struct cmt_mgmt *cm, TPnode *tpnode) {
    uint32_t pos = tp_hash(tpnode->tvpn);
    QTAILQ_INSERT_HEAD(&cm->hash_tp_table[pos], tpnode, h_entry);
}

static void remove_tpnode(struct cmt_mgmt *cm, TPnode *tpnode) {
    ftl_assert(tpnode->cmt_entry_cnt == 0 && lm_cnt == 0);
    QTAILQ_REMOVE(&cm->TPnode_list, tpnode, lru_entry); /* detach from the lru list*/
    uint32_t pos = tp_hash(tpnode->tvpn);
    QTAILQ_REMOVE(&cm->hash_tp_table[pos], tpnode, h_entry); /* detach from the hash table */
    cm->live_tpnode_cnt--;
    cm->free_cmt_entry_cnt++;
}

/*
* Only reclaim cmt_entry from the cold list.
*/
static void reclaim_cmt_entry(struct cmt_mgmt *cm, TPnode *tpnode, cmt_entry *mapping_entry) {
    QTAILQ_REMOVE(&tpnode->cold_list, mapping_entry, entry); /* detach from the lru list*/
    tpnode->cmt_entry_cnt--;
    uint64_t pos = cmt_hash(mapping_entry->single.lpn);
    QTAILQ_REMOVE(&cm->hash_tp_table[pos], mapping_entry, h_entry); /* detach from the hash table */
    cm->used_cmt_entry_cnt--;
    cm->free_cmt_entry_cnt++;
    QTAILQ_INSERT_HEAD(&cm->free_cmt_entry_list, mapping_entry, entry);
    if (tpnode->cmt_entry_cnt == 0 && tpnode->lm_cnt == 0) {
        remove_tpnode(cm, tpnode);
    }
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
            ppn_offset[traindata_idx] = mapping_entry->single.ppn.g.lg_offset;
            traindata_idx++;
        }
    }
    ftl_assert(traindata_idx == num - 1);
    QsortKV(lpn_offset, ppn_offset, 0, num);

    // fit the model
    struct ransac_t model = {.bestfit_cnt = 0, .maxerror = 1.0f, .samplesize = num/10};
    RANSAC_fit(lpn_offset, ppn_offset, num, &model);

    // mark the precise points in gtd-granularity
    // get the model range [start, start+len)
    linear_model_t *lm = &model.bestmodel;
    uint32_t lm_start = 0, lm_len = 0; 
    for (uint32_t i = 0; i < num; i++) {
        if (roundf(lm->w * lpn_offset[i] + lm->b) == ppn_offset[i]) {
            precise_bitmap[i] = true;          
        } else {
            precise_bitmap[i] = false;
        };
    }
    for (uint32_t i = 0; i < num; i++) {
        if ()
    }
}