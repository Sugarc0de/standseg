/*
 * linux_compat.c
 *
 * main.c computes
 *
 *   sproc.reclaim_trigger = (getpagesize() * MIN_RECLAIM_PAGES) /
 *       (sizeof(region) + sizeof(neighbor) + sizeof(float) * sproc.nbands) + 1;
 *
 * i.e. the number of dead regions that must accumulate before
 * compact_region_list() runs.  Compaction RENUMBERS every region, and the
 * merge loop walks regions in id order, so the trigger value is part of the
 * algorithm's output, not a performance knob.
 *
 * getpagesize() is 4096 on the Linux/x86-64 machine that produced the golden
 * fixtures but 16384 on macOS/arm64, which would give 3641 instead of 911 and
 * change which passes compact.  (Derived independently from
 * tests/golden/test_3456/expected/myseg.log: the observed compaction pattern
 * pins the trigger to the interval (869, 939]; 4096*8/36+1 = 911.)
 *
 * Defining getpagesize() here makes our definition win over libc's at link
 * time, exactly as glibc_random.c does for random().
 */

int
getpagesize(void)
{
    return 4096;
}
