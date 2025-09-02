#include <htslib/sam.h>
#include "bwa.h"
#include "stdio.h"
#include <assert.h>
#include <string.h>
#include <stdbool.h>

bwaidx_t *idx;
char *filename1;

char *filename2;

static inline char nt16_to_ACGTN(uint8_t x)
{
    switch (x)
    {
    case 1:
        return 'A';
    case 2:
        return 'C';
    case 4:
        return 'G';
    case 8:
        return 'T';
    case 15:
        return 'N';
    default:
        printf("Invalid nt16 code: %d\n", x);
        abort();
        return 'N';
    }
}
static inline char comp_base(char b)
{
    switch (b)
    {
    case 'A':
        return 'T';
    case 'C':
        return 'G';
    case 'G':
        return 'C';
    case 'T':
        return 'A';
    default:
        return 'N';
    }
}

int get_ref_len(uint32_t n_cigar, uint32_t *cigar,int *ref_len, int *read_len){
    *ref_len=0;
    *read_len=0;
    for(int i=0;i<n_cigar;i++){
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);
        switch (op){
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                *ref_len+=len;
                *read_len+=len;
                break;
            case BAM_CINS:
                *read_len+=len;
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                *ref_len+=len;
                break;
            case BAM_CSOFT_CLIP:
                *read_len+=len;
                break;
            case BAM_CHARD_CLIP:
                break;
            case BAM_CPAD:
                break;
            default:
                fprintf(stderr,"Invalid CIGAR op: %u\n",op);
                return -1;
        }
    }
    return 0;
}

int cigar_seq_count(uint64_t *cigar_count,uint32_t n_cigar, uint32_t *cigar){
    for(int i=0;i<n_cigar;i++){
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);
        cigar_count[op]+=len;
        switch (op){
            case BAM_CMATCH://BAM_CMATCH 在 BWA 中表示匹配或不匹配，需要特别处理
                // cigar_count[BAM_CMATCH]+=len;
                break;
            //BWA中没有使用 BAM_CEQUAL 和 BAM_CDIFF，而是全部用 BAM_CMATCH
            // case BAM_CEQUAL:
            //     cigar_count[BAM_CEQUAL]+=len;
            //     break;
            // case BAM_CDIFF:
            //     cigar_count[BAM_CDIFF]+=len;
            //     break;
            case BAM_CINS:
                cigar_count[BAM_CINS]+=len;
                break;
            case BAM_CDEL:
                cigar_count[BAM_CDEL]+=len;
                break;
            case BAM_CREF_SKIP:
                cigar_count[BAM_CREF_SKIP]+=len;
                break;
            case BAM_CSOFT_CLIP:
                cigar_count[BAM_CSOFT_CLIP]+=len;
                break;
            case BAM_CHARD_CLIP:
                cigar_count[BAM_CHARD_CLIP]+=len;
                break;
            case BAM_CPAD:
                cigar_count[BAM_CPAD]+=len;
                break;
            default:
                fprintf(stderr,"Invalid CIGAR op: %u\n",op);
                return -1;
        }
    }
    return 0;
}

//计算 CIGAR 是否与参考序列匹配 ，统计匹配和不匹配的碱基数
int check_cigar_count(uint64_t *cigar_count,uint32_t n_cigar, uint32_t *cigar, int64_t pos /*0-based in contig*/,
    uint8_t *read_seq, int seq_len, int is_rev, int tid, bam_hdr_t *header)
{
    int read_pos = 0;
    uint8_t *pac = idx->pac;
    uint8_t *rseq;
    int64_t rlen;

    // compute global pac position: add contig offset in pac
    int64_t contig_offset = idx->bns->anns[tid].offset; // 注意字段名与 BWA 版本一致性
    int64_t gpos = contig_offset + pos;                 // gpos 是 pac 的 0-based 全局坐标

    for (int i = 0; i < n_cigar; i++)
    {
        uint32_t op = bam_cigar_op(cigar[i]);
        uint32_t len = bam_cigar_oplen(cigar[i]);

        switch (op)
        {
        case BAM_CMATCH:
            // bns_get_seq 的 end 参数是否闭区间请根据你 BWA 版本确认；此处用 [gpos, gpos+len)
            rseq = bns_get_seq(idx->bns->l_pac, pac, gpos, gpos + len, &rlen);
            if (rlen != len)
            {
                fprintf(stderr, "DEBUG: rlen(%ld) != len(%u) for tid=%d (%s) gpos=%ld\n",
                        rlen, len, tid, header->target_name[tid], gpos);
                free(rseq);
                return -1;
            }
            
            for (int j = 0; j < len; j++)
            {
                int ref_code = (int)rseq[j];
                char ref_base = "ACGTN"[ref_code];
                int qidx = is_rev ? (seq_len - 1 - (read_pos + j)) : (read_pos + j);
                uint8_t nt16 = bam_seqi(read_seq, qidx);
                char read_base = nt16_to_ACGTN(nt16);
                // if (is_rev)i
                //     read_base = comp_base(read_base);
                if (read_base == 'N') cigar_count[BAM_CMATCH]++; 
                if (ref_base != read_base)
                {
                    cigar_count[BAM_CDIFF]++;
                }else{
                    cigar_count[BAM_CEQUAL]++;
                }
            }
            free(rseq);
            gpos += len; // advance global position
            pos += len;  // and contig-local pos
            read_pos += len;
            break;
        case BAM_CEQUAL:
        case BAM_CDIFF:
            break;
        case BAM_CINS:
            read_pos += len;
            break;
        case BAM_CDEL:
        case BAM_CREF_SKIP:
            gpos += len;
            pos += len;
            break;
        case BAM_CSOFT_CLIP:
            read_pos += len;
            break;
        case BAM_CHARD_CLIP:
            break;
        case BAM_CPAD:
            break;
        default:
            return -1;
        }
    }
    return 0;
}


// 替换用的核心函数：直接从 idx->pac 提取参考碱基（2-bit packed）
static inline char pac_code_to_base(uint8_t code){
    // pac 存的是 0->A,1->C,2->G,3->T (你以前用的顺序是 ACGT)
    // 如果你的 pac 定义不同请调整映射
    switch(code & 3){
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

void count_cigar_in_file(samFile *fp,bam_hdr_t *header,bam1_t *aln,uint64_t *cigar_count,uint64_t *map_seq_count,int64_t max_read_num){
    int read_num=0;
    while (1)
    {

        if (sam_read1(fp, header, aln) < 0)
            break;
        if(max_read_num>0){
            read_num++;
            if(read_num>max_read_num)break; 
        }
        const bam1_core_t *c = &aln->core;
        int is_rev = c->flag & BAM_FREVERSE;
        int tid = c->tid;
        if(c->flag&BAM_FUNMAP){
            map_seq_count[0]++;
            continue;
        }
            
        if (tid < 0){
            map_seq_count[0]++;
            continue;
        }
            
        if (cigar_seq_count(cigar_count,c->n_cigar, bam_get_cigar(aln)) < 0)
        {
            printf("file cigar count failed for read %s:\n", bam_get_qname(aln));
            printf("flag: %d,pos: %ld ", c->flag, c->pos);
            return;
        }
        if(check_cigar_count(cigar_count,c->n_cigar, bam_get_cigar(aln), c->pos, bam_get_seq(aln), c->l_qseq, is_rev, tid, header) < 0){
            printf("file cigar check failed for read %s:\n", bam_get_qname(aln));
            printf("flag: %d,pos: %ld ", c->flag, c->pos);
            return;
        }
        map_seq_count[1]++;
    }
}

void test(samFile *fp1,bam_hdr_t *header1,bam1_t *aln1,samFile *fp2,bam_hdr_t *header2,bam1_t *aln2){
    uint64_t cigar_count_f1[9] = {0};
    uint64_t cigar_count_f2[9] = {0};
    uint64_t map_seq_count_f1[2] = {0};//0 for unmapped, 1 for mapped
    uint64_t map_seq_count_f2[2] = {0};//0 for unmapped, 1 for mapped
    //在这里对CIGAR统计上，由于BWA中没有使用 BAM_CEQUAL 和 BAM_CDIFF，而是全部用 BAM_CMATCH，因此计数使用中BAM_CEQUAL 和 BAM_CDIFF位置统计 BAM_CMATCH比对上和每比对上的数量，原BAM_CMATCH位置统计碱基为N的数量
    count_cigar_in_file(fp1,header1,aln1,cigar_count_f1,map_seq_count_f1,0);
    count_cigar_in_file(fp2,header2,aln2,cigar_count_f2,map_seq_count_f2,0);
    printf("In this CIGAR statistics, since BAM_CEQUAL and BAM_CDIFF were not used in BWA but instead all were counted as BAM_CMATCH, the count used for BAM_CEQUAL and BAM_CDIFF is based on the statistics of BAM_CMATCH alignments and the number of each alignment, while the original BAM_CMATCH position counts the number of bases as N.\n");
    printf("CIGAR op seq counts for file1:\n");
    printf("M: %lu,%lu\n",cigar_count_f1[BAM_CMATCH],cigar_count_f2[BAM_CMATCH]);
    printf("=: %lu,%lu\n",cigar_count_f1[BAM_CEQUAL],cigar_count_f2[BAM_CEQUAL]);
    printf("X: %lu,%lu\n",cigar_count_f1[BAM_CDIFF],cigar_count_f2[BAM_CDIFF]);
    printf("I: %lu,%lu\n",cigar_count_f1[BAM_CINS],cigar_count_f2[BAM_CINS]);
    printf("D: %lu,%lu\n",cigar_count_f1[BAM_CDEL],cigar_count_f2[BAM_CDEL]);
    printf("N: %lu,%lu\n",cigar_count_f1[BAM_CREF_SKIP],cigar_count_f2[BAM_CREF_SKIP]);
    printf("S: %lu,%lu\n",cigar_count_f1[BAM_CSOFT_CLIP],cigar_count_f2[BAM_CSOFT_CLIP]);
    printf("H: %lu,%lu\n",cigar_count_f1[BAM_CHARD_CLIP],cigar_count_f2[BAM_CHARD_CLIP]);
    printf("P: %lu,%lu\n",cigar_count_f1[BAM_CPAD],cigar_count_f2[BAM_CPAD]);
    printf("Mapped and unmapped reads:\n");
    printf("File1: mapped: %lu, unmapped: %lu\n",map_seq_count_f1[1],map_seq_count_f1[0]);
    printf("File2: mapped: %lu, unmapped: %lu\n",map_seq_count_f2[1],map_seq_count_f2[0]);
}



int main(int argc ,char**argv){
    if(argc<4){
        printf("Usage: ./a.out <bwa_index> <sam_file1> <sam_file2>\n");
        return 1;
    }
    idx = bwa_idx_load(argv[1], BWA_IDX_ALL);
    filename1 = argv[2];
    filename2 = argv[3];
    samFile *fp1 = sam_open(filename1, "r");
    if (fp1 == NULL) {
        perror("无法打开文件");
        return 1;
    }
    samFile *fp2 = sam_open(filename2, "r");
    if (fp2 == NULL) {
        perror("无法打开文件");
        return 1;
    }
    bam_hdr_t *header1 = sam_hdr_read(fp1);
    bam_hdr_t *header2 = sam_hdr_read(fp2);
    bam1_t *aln1 = bam_init1();
    bam1_t *aln2 = bam_init1();
    printf("NEXT TO TEST\n");
    test(fp1,header1,aln1,fp2,header2,aln2);
    bam_destroy1(aln1);
    bam_destroy1(aln2);
    bam_hdr_destroy(header1);
    bam_hdr_destroy(header2);
    sam_close(fp1);
    sam_close(fp2);
    bwa_idx_destroy(idx);
    return 0;
}