#ifndef __FEMU_TPFTL_H
#define __FEMU_TPFTL_H

#include "../nvme.h"

#define INVALID_PPA     (~(0ULL))
#define INVALID_LPN     (~(0ULL))
#define UNMAPPED_PPA    (~(0ULL))

#define CMT_HASH_SIZE (24593ULL)
#define TP_HASH_SIZE (24593ULL)
#define LG_HASH_SIZE (24593ULL)
#define SUB_SPACE_SIZE (1ULL << 28) /* size in term of bytes, 256MB */

#define ENT_PER_TP (512ULL)

enum {
    NAND_READ =  0,
    NAND_WRITE = 1,
    NAND_ERASE = 2,

    NAND_READ_LATENCY = 40000,
    NAND_PROG_LATENCY = 200000,
    NAND_ERASE_LATENCY = 2000000,
};

enum {
    USER_IO = 0,
    GC_IO = 1,
};

enum {
    SEC_FREE = 0,
    SEC_INVALID = 1,
    SEC_VALID = 2,

    PG_FREE = 0,
    PG_INVALID = 1,
    PG_VALID = 2
};

enum {
    CLEAN = 0,
    DIRTY = 1
};

enum {
    HEAD = 0,
    TAIL = 1
};

enum {
    NONE = 0,
    DATA = 1,
    TRANS = 2
};

enum {
    SINGLE_MAPPING = 0,
    LINEAR_MODEL = 1
};

enum {
    FEMU_ENABLE_GC_DELAY = 1,
    FEMU_DISABLE_GC_DELAY = 2,

    FEMU_ENABLE_DELAY_EMU = 3,
    FEMU_DISABLE_DELAY_EMU = 4,

    FEMU_RESET_ACCT = 5,
    FEMU_ENABLE_LOG = 6,
    FEMU_DISABLE_LOG = 7,

    FEMU_RESET_STAT = 8,
    FEMU_PRINT_STAT = 9,
};


#define BLK_BITS    (16)
#define PG_BITS     (16)
#define SEC_BITS    (8)
#define PL_BITS     (8)
#define LUN_BITS    (8)
#define CH_BITS     (7)

/* describe a physical page addr */
struct ppa {
    union {
        struct {
            uint64_t blk : BLK_BITS;
            uint64_t pg  : PG_BITS;
            uint64_t sec : SEC_BITS;
            uint64_t pl  : PL_BITS;
            uint64_t lun : LUN_BITS;
            uint64_t ch  : CH_BITS;
            uint64_t rsv : 1;
        } g;

        uint64_t ppa;
    };
};

typedef int nand_sec_status_t;

struct nand_page {
    nand_sec_status_t *sec;
    int nsecs;
    int status;
};

struct nand_block {
    struct nand_page *pg;
    struct ppa blk_ppa;
    bool is_allocated;
    bool is_bad;
    int npgs;
    int ipc; /* invalid page count */
    int vpc; /* valid page count */
    int erase_cnt;
    QTAILQ_ENTRY(nand_block) entry; /* in free block list */
};

struct nand_plane {
    struct nand_block *blk;
    int nblks;
    int nfreeblks;
    QTAILQ_HEAD(free_blk_list, nand_block) free_blk_list;
};

struct nand_lun {
    struct nand_plane *pl;
    int npls;
    uint64_t next_lun_avail_time;
    bool busy;
    uint64_t gc_endtime;
};

struct ssd_channel {
    struct nand_lun *lun;
    int nluns;
    uint64_t next_ch_avail_time;
    bool busy;
    uint64_t gc_endtime;
};

struct ssdparams {
    int secsz;        /* sector size in bytes */
    int secs_per_pg;  /* # of sectors per page */
    int pgs_per_blk;  /* # of NAND pages per block */
    int blks_per_pl;  /* # of blocks per plane */
    int pls_per_lun;  /* # of planes per LUN (Die) */
    int luns_per_ch;  /* # of LUNs per channel */
    int nchs;         /* # of channels in the SSD */

    int pg_rd_lat;    /* NAND page read latency in nanoseconds */
    int pg_wr_lat;    /* NAND page program latency in nanoseconds */
    int blk_er_lat;   /* NAND block erase latency in nanoseconds */
    int ch_xfer_lat;  /* channel transfer latency for one page in nanoseconds
                       * this defines the channel bandwith
                       */

    double gc_thres_pcent;
    int gc_thres_lines;
    double gc_thres_pcent_high;
    int gc_thres_lines_high;
    bool enable_gc_delay;

    /* below are all calculated values */
    int secs_per_blk; /* # of sectors per block */
    int secs_per_pl;  /* # of sectors per plane */
    int secs_per_lun; /* # of sectors per LUN */
    int secs_per_ch;  /* # of sectors per channel */
    int tt_secs;      /* # of sectors in the SSD */

    int pgs_per_pl;   /* # of pages per plane */
    int pgs_per_lun;  /* # of pages per LUN (Die) */
    int pgs_per_ch;   /* # of pages per channel */
    int tt_pgs;       /* total # of pages in the SSD */

    int blks_per_lun; /* # of blocks per LUN */
    int blks_per_ch;  /* # of blocks per channel */
    int tt_blks;      /* total # of blocks in the SSD */

    int blks_per_lg;  /* length of linear group */

    int pls_per_ch;   /* # of planes per channel */
    int tt_pls;       /* total # of planes in the SSD */

    int tt_luns;      /* total # of LUNs in the SSD */

    int tt_subspaces;  /* total # of sub spaces in the SSD */

    int ents_per_pg;
    int tt_cmt_size;
    int tt_gtd_size;
    bool enable_request_prefetch;
    bool enable_select_prefetch;
};

struct nand_cmd {
    int type;
    int cmd;
    int64_t stime; /* Coperd: request arrival time */
};

/**
 * @brief cmt entry struct and cmt TPnode struct and cmt management struct
 * 
 */
typedef struct cmt_entry {
    int type;
    int dirty;
    bool prefetch;
    union {
        struct {
            uint64_t lpn;
            uint64_t ppn;
        } single;
        struct {
            float slope;
            float intercept;
            uint32_t start_lpn_offset;
            uint32_t lg_id;
        } lm;
    };
    QTAILQ_ENTRY(cmt_entry) entry; /* for per-TPnode list or free list*/
    QTAILQ_ENTRY(cmt_entry) h_entry; /* for hash */
    uint64_t next_avail_time; 
} cmt_entry;

typedef struct TPnode {
    uint64_t tvpn;
    int cmt_entry_cnt;
    int lm_cnt;
    QTAILQ_ENTRY(TPnode) lru_entry;
    QTAILQ_ENTRY(TPnode) h_entry;
    QTAILQ_HEAD(lm_list, cmt_entry) lm_list; /* cmt entry - head insert */
    QTAILQ_HEAD(cmt_entry_list, cmt_entry) cmt_entry_list; /* cmt entry - head insert */
} TPnode;

/*
* The cached mapping table is to serve unlinaer mappings
*/
struct cmt_mgmt {
    cmt_entry *cmt_entries;
    QTAILQ_HEAD(free_cmt_entry_list, cmt_entry) free_cmt_entry_list;
    //QTAIQ_HEAD 为热度最高的
    QTAILQ_HEAD(TPnode_list, TPnode) TPnode_list;
    int tt_TPnodes;
    int tt_entries;
    int free_cmt_entry_cnt;
    int used_cmt_entry_cnt;
    int live_tpnode_cnt;
    // use for selective prefetching;
    int counter;

    QTAILQ_HEAD(cmt_hash_table, cmt_entry) hash_mapping_table[CMT_HASH_SIZE];
    QTAILQ_HEAD(tp_table, TPnode) hash_tp_table[TP_HASH_SIZE]; /* virtual translation page number -> TPnode*/
};

struct lg_write_pointer {
    int pg;
    int blk;
};

struct lg_allocate_pointer {
    int ch;
    int lun;
    int pl;
};

/*
* [action] open -> close: liner-group mapping table -> linear model
* group inner offset -> real blk id
* create liner group: allocate blocks
*/
struct linear_group {
    bool is_open;
    int id; // the unique linear group id
    int len; // number of blks
    int ipc; /* invalid page count in this line */
    int vpc; /* valid page count in this line */
    /* position in the priority queue for victim lines */
    size_t pos;

    int type;
    struct ppa *blks; // the vector recording blks' addresses
    uint16_t *reverse_lpns; /* TODO: discard the additional reverse lpn records */
    uint32_t start_offset; /* the start offset of un trained records */
    struct lg_write_pointer lg_wp; // per linear group write pointer

    QTAILQ_ENTRY(linear_group) entry; // for the chain of the same sub space
};

/*
* [map] lpn -> linear groups (the first one is open)
* check lun free blk list 
*/
struct lg_mgmt {
    QTAILQ_HEAD(lg_list, linear_group) lg_list[LG_HASH_SIZE]; // hash(id) -> linear groups. !! the additional one is for translation pages
    pqueue_t *victim_line_pq; // closed linear group having invalid page becomes viticm member
    struct lg_allocate_pointer lg_ap;
    int tt_lg;
    int open_lg_cnt;
    int close_lg_cnt;
    int victim_lg_cnt;
};

struct statistics {
    uint64_t cmt_hit_cnt;
    uint64_t cmt_miss_cnt;
    double cmt_hit_ratio;
    uint64_t access_cnt;
};

struct ssd {
    char *ssdname;
    struct ssdparams sp;
    struct ssd_channel *ch;
    struct ppa *maptbl; /* page level mapping table */
    uint64_t *rmap;     /* reverse mapptbl, assume it's stored in OOB */
    struct lg_mgmt lgm;
    
    struct cmt_mgmt cm;
    struct ppa *gtd;

    /* lockless ring for communication with NVMe IO thread */
    struct rte_ring **to_ftl;
    struct rte_ring **to_poller;
    bool *dataplane_started_ptr;
    QemuThread ftl_thread;

    struct statistics stat;
    FILE *fpr, *fpw;
};

void ssd_init(FemuCtrl *n);

#ifdef FEMU_DEBUG_FTL
#define ftl_debug(fmt, ...) \
    do { printf("[FEMU] FTL-Dbg: " fmt, ## __VA_ARGS__); } while (0)
#else
#define ftl_debug(fmt, ...) \
    do { } while (0)
#endif

#define ftl_err(fmt, ...) \
    do { fprintf(stderr, "[FEMU] FTL-Err: " fmt, ## __VA_ARGS__); } while (0)

#define ftl_log(fmt, ...) \
    do { printf("[FEMU] FTL-Log: " fmt, ## __VA_ARGS__); } while (0)


/* FEMU assert() */
#ifdef FEMU_DEBUG_FTL
#define ftl_assert(expression) assert(expression)
#else
#define ftl_assert(expression)
#endif

#endif
