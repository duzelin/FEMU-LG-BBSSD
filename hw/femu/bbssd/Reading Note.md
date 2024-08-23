# Reading Note
Author: Du Zelin

## FTL

enum NAND latency
enum user/gc io
enum sector/page state
enum FEMU related settings
struct ppa
channel -> lun -> plane -> block -> page -> sectors
struct ssdparams
struct write pointer
struct line(block) mgmt
ssd: ssdparams, channels, mapping table, reverse mapping table(OOB), write pointer, line manager, NVMe IO ring pair, QemuThread ftl
void function ssd_init()

ftl_thread
    - dequeue ftl ring
    - process cmd
        - ssd_write()
        - ssd_read()
    - enqueue poller ring
    - garbage collection
## DFTL
+ struct cmt mgmt
+ struct trans_write_pointer
+ process_translation_page_write()
## LeaFTL

## Learned FTL
