#include "lgftl.h"

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
        inliers += (error[i] < maxerror ? 1 : 0);
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
    
    g_free(sample_x);
    g_free(sample_y);
    g_free(error);
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

static uint64_t virtual_tpage_num (struct ppa *gtd, uint64_t lpn) {
    uint64_t gtd_offset = lpn / ENT_PER_TP;
    while (gtd_offset && gtd[gtd_offset].ppa == gtd[gtd_offset - 1].ppa) {
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

static struct TPnode *find_tpnode(struct cmt_mgmt *cm, struct ppa *gtd, uint64_t lpn) {
    uint64_t tvpn = virtual_tpage_num(gtd, lpn);
    uint32_t pos = tp_hash(tvpn);
    TPnode *tpnode = NULL;
    QTAILQ_FOREACH(tpnode, &cm->hash_tp_table[pos], h_entry) {
        if (tpnode->tvpn == tvpn) {
            break;
        }
    }
    return tpnode;
}

static struct cmt_entry *find_lm_entry(struct cmt_mgmt *cm, struct ppa *gtd, uint64_t lpn) {
    uint32_t offset_in_subspace = lpn % SUB_SPACE_SIZE;
    struct TPnode *tpnode = find_tpnode(cm, gtd, lpn);
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

static void reclaim_cmt_entry(struct cmt_mgmt *cm, TPnode *tpnode, cmt_entry *mapping_entry) {
    QTAILQ_REMOVE(&tpnode->cmt_entry_list, mapping_entry, entry); /* detach from the lru list*/
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
    new_lg->reverse_lpns = g_malloc0(sizeof(uint16_t) * new_lg->len);

    QTAILQ_INSERT_HEAD(&lgm->lg_list[sub_space_id], new_lg, entry);
    lgm->tt_lg++;
    lgm->open_lg_cnt++;
}

void close_lg(struct ssd *ssd, struct linear_group *lg) {
    uint32_t start_offset = lg->start_offset;
    uint32_t final_offset = lg->lg_wp.pg * lg->len + lg->lg_wp.blk;
    for (int i = start_offset; i < final_offset; i++) {
        
    }
}