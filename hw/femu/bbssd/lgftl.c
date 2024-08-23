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

static bool create_linear_group(struct ssd *ssd) {
    struct ssd_channel* channels = ssd->ch;

    for (int i = 0; i < channels->nluns; i++) {
        
    }
}