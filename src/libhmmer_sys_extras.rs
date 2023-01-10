#![allow(non_upper_case_globals)]

// /* Return codes for error handler
//  */
// /*::cexcerpt::statuscodes::begin::*/
// #define eslOK              0    /* no error/success             */
// #define eslFAIL            1    /* failure                      */
// #define eslEOL             2    /* end-of-line (often normal)   */
// #define eslEOF             3    /* end-of-file (often normal)   */
// #define eslEOD             4    /* end-of-data (often normal)   */
// #define eslEMEM            5    /* malloc or realloc failed     */
// #define eslENOTFOUND       6    /* file or key not found        */
// #define eslEFORMAT         7    /* file format not correct      */
// #define eslEAMBIGUOUS      8    /* an ambiguity of some sort    */
// #define eslEDIVZERO        9    /* attempted div by zero        */
// #define eslEINCOMPAT      10    /* incompatible parameters      */
// #define eslEINVAL         11    /* invalid argument/parameter   */
// #define eslESYS           12    /* generic system call failure  */
// #define eslECORRUPT       13    /* unexpected data corruption   */
// #define eslEINCONCEIVABLE 14    /* "can't happen" error         */
// #define eslESYNTAX        15    /* invalid user input syntax    */
// #define eslERANGE         16    /* value out of allowed range   */
// #define eslEDUP           17    /* saw a duplicate of something */
// #define eslENOHALT        18    /* a failure to converge        */      
// #define eslENORESULT      19    /* no result was obtained       */
// #define eslENODATA        20    /* no data provided, file empty */
// #define eslETYPE          21    /* invalid type of argument     */
// #define eslEOVERWRITE     22    /* attempted to overwrite data  */
// #define eslENOSPACE       23    /* ran out of some resource     */
// #define eslEUNIMPLEMENTED 24    /* feature is unimplemented     */
// #define eslENOFORMAT      25	/* couldn't guess file format   */
// #define eslENOALPHABET    26	/* couldn't guess seq alphabet  */
// #define eslEWRITE         27   	/* write failed (fprintf, etc)  */
// #define eslEINACCURATE    28    /* return val may be inaccurate */
// /*::cexcerpt::statuscodes::end::*/
pub const eslOK: i32 = 0;
pub const eslFAIL: i32 = 1;
pub const eslEOL: i32 = 2;
pub const eslEOF: i32 = 3;
pub const eslEOD: i32 = 4;
pub const eslEMEM: i32 = 5;
pub const eslENOTFOUND: i32 = 6;
pub const eslEFORMAT: i32 = 7;
pub const eslEAMBIGUOUS: i32 = 8;
pub const eslEDIVZERO: i32 = 9;
pub const eslEINCOMPAT: i32 = 10;
pub const eslEINVAL: i32 = 11;
pub const eslESYS: i32 = 12;
pub const eslECORRUPT: i32 = 13;
pub const eslEINCONCEIVABLE: i32 = 14;
pub const eslESYNTAX: i32 = 15;
pub const eslERANGE: i32 = 16;
pub const eslEDUP: i32 = 17;
pub const eslENOHALT: i32 = 18;
pub const eslENORESULT: i32 = 19;
pub const eslENODATA: i32 = 20;
pub const eslETYPE: i32 = 21;
pub const eslEOVERWRITE: i32 = 22;
pub const eslENOSPACE: i32 = 23;
pub const eslEUNIMPLEMENTED: i32 = 24;
pub const eslENOFORMAT: i32 = 25;
pub const eslENOALPHABET: i32 = 26;
pub const eslEWRITE: i32 = 27;
pub const eslEINACCURATE: i32 = 28;