#![allow(non_upper_case_globals)]

// /*****************************************************************
//  * 16. P7_PIPELINE: H3's accelerated seq/profile comparison pipeline
//  *****************************************************************/
//  enum p7_pipemodes_e { p7_SEARCH_SEQS = 0, p7_SCAN_MODELS = 1 };
//  enum p7_zsetby_e    { p7_ZSETBY_NTARGETS = 0, p7_ZSETBY_OPTION = 1, p7_ZSETBY_FILEINFO = 2 };
//  enum p7_complementarity_e { p7_NOCOMPLEMENT    = 0, p7_COMPLEMENT   = 1 };
pub const p7_SEARCH_SEQS: i32 = 0;
// TODO: Add more, but eh for now.

// #define p7_LOCAL     1		/* multihit local:  "fs" mode   */
pub const p7_LOCAL: i32 = 1;
