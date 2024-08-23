#include "lgftl.h"

static void *ftl_thread(void *arg);

/* process hash */
static inline uint64_t cmt_hash(uint64_t lpn)
{
    return lpn % CMT_HASH_SIZE;
}

static inline uint64_t tp_hash(uint64_t tvpn)
{
    return tvpn % TP_HASH_SIZE;
}

static struct cmt_entry* find_hash_entry(hash_table *ht, uint64_t lpn)
{
    uint64_t pos = cmt_hash(lpn);
    cmt_entry *entry = ht->cmt_table[pos];
    while (entry != NULL && entry->lpn != lpn) {
        entry = entry->next;
    }
    return entry;
}

static struct TPnode* find_hash_tpnode(hash_table *ht, uint64_t tvpn)
{
    uint64_t pos = tp_hash(tvpn);
    TPnode *tpnode = ht->tp_table[pos];
    while (tpnode != NULL && tpnode->tvpn != tvpn) {
        tpnode = tpnode->next;
    }
    return tpnode;
}

static void insert_cmt_hashtable(hash_table *ht, cmt_entry *entry) 
{
    uint64_t pos = cmt_hash(entry->lpn);
    entry->next = ht->cmt_table[pos];
    ht->cmt_table[pos] = entry;
}

static void insert_tp_hashtable(hash_table *ht, TPnode *tpnode) 
{
    uint64_t pos = tp_hash(tpnode->tvpn);
    tpnode->next = ht->tp_table[pos];
    ht->tp_table[pos] = tpnode;
}

static bool delete_cmt_hashnode(hash_table *ht, cmt_entry *entry)
{
    uint64_t pos = cmt_hash(entry->lpn);
    cmt_entry *tmp_entry = ht->cmt_table[pos], *pre_entry;
    if (tmp_entry == entry) {
        ht->cmt_table[pos] = tmp_entry->next;
        tmp_entry->next = NULL;
    } else {
        pre_entry = tmp_entry;
        tmp_entry = tmp_entry->next;
        while (tmp_entry != NULL && tmp_entry != entry) {
            pre_entry = tmp_entry;
            tmp_entry = tmp_entry->next;
        }
        if (tmp_entry == NULL)
            return false;
        pre_entry->next = tmp_entry->next;
        tmp_entry->next = NULL;
    }
    return true;
}

static bool delete_tp_hashnode(hash_table *ht, TPnode *tpnode)
{
    uint64_t pos = tp_hash(tpnode->tvpn);
    TPnode *tmp_tp = ht->tp_table[pos], *pre_tp;
    // The single tpnode in the target hash slot.
    if (tmp_tp == tpnode) {
        ht->tp_table[pos] = tmp_tp->next;
        tmp_tp->next = NULL;
    } else {
        pre_tp = tmp_tp;
        tmp_tp = tmp_tp->next;
        while (tmp_tp != NULL && tmp_tp != tpnode) {
            pre_tp = tmp_tp;
            tmp_tp = tmp_tp->next;
        }
        if (tmp_tp == NULL)
            return false;
        pre_tp->next = tmp_tp->next;
        tmp_tp->next = NULL;
    }
    return true;
}

struct ppa allocate_blk(struct ssd *ssd) {
    struct lg_allocate_pointer* lg_ap = ssd->lgm.lg_ap;
    struct ssd_channel* channels = ssd->ch;
    int ch_num = lg_ap->ch;
    int lun_num = lg_ap->lun;
    int pl_num = lg_ap->pl;
    struct nand_plane* target_plane = &channels[ch_num].lun[lun_num].pl[pl_num];
    struct nand_block* candidate_block = NULL;

    candidate_block = QTAILQ_FIRST(&target_plane->free_blk_list);
        if (!candidate_block) {
        return NULL;
    }

    QTAILQ_REMOVE(&target_plane->free_blk_list, candidate_block, entry);  
    candidate_block.is_allocated = true;   
    target_plane->nfreeblks--;
    return candidate_block->blk_ppa;
}

static void advance_lg_allocate_pointer(struct ssd *ssd) {
    struct lg_allocate_pointer* lg_ap = ssd->lgm.lg_ap;
    struct ssdparams* sp = ssd->sp;
    lg_ap.ch++;
    if(lg_ap.ch == sp->nchs) {
        lg_ap.ch = 0;
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

struct ppa* allocate_blks_for_lg(struct ssd *ssd) {
    struct ssdparams* sp = ssd->sp;
    struct ssd_channel* channels = ssd->ch;
    struct lg_mgmt* lgm= ssd->lgm;

    int lg_len = sp->blks_per_lg;
    struct ppa* blks = g_malloc0(sizeof(struct ppa) * lg_len);
    for (int i = 0; i < lg_len; i++) {
        // add one blk from current lun to new linear group
        blks[i] = allocate_blk(ssd);
        advance_lg_allocate_pointer(ssd);
        
    }
}