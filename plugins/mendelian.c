/* The MIT License

   Copyright (c) 2015-2022 Genome Research Ltd.

   Author: Petr Danecek <pd3@sanger.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <getopt.h>
#include <math.h>
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>     // for isatty
#include "../bcftools.h"
#include "../regidx.h"

#define MODE_COUNT     1
#define MODE_LIST_GOOD 2
#define MODE_LIST_BAD  4
#define MODE_DELETE    8
#define MODE_ANNOTATE  16
#define MODE_LIST_SKIP 32

typedef struct
{
    int nok, nbad;
    int imother,ifather,ichild;
}
trio_t;

typedef struct
{
    int mpl, fpl, cpl;  // ploidies - mother, father, child
    int mal, fal;       // expect an allele from mother and father
}
rule_t;

typedef struct _args_t
{
    regidx_t *rules;
    regitr_t *itr, *itr_ori;
    bcf_hdr_t *hdr;
    htsFile *out_fh;
    int32_t *gt_arr;
    int mode;
    int ngt_arr, nrec;
    trio_t *trios;
    int ntrios, mtrios;
    int output_type, clevel;
    char *output_fname;
    bcf_srs_t *sr;
}
args_t;

static args_t args;
static int parse_rules(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr);
static bcf1_t *process(bcf1_t *rec);

const char *about(void)
{
    return "Count Mendelian consistent / inconsistent genotypes [DEPRECATED, use mendelian2 instead]\n";
}

typedef struct
{
    const char *alias, *about, *rules;
}
rules_predef_t;

static rules_predef_t rules_predefs[] =
{
    { .alias = "GRCh37",
      .about = "Human Genome reference assembly GRCh37 / hg19, both chr naming conventions",
      .rules =
            "   X:1-60000               M/M + F > M\n"
            "   X:1-60000               M/M + F > M/F\n"
            "   X:2699521-154931043     M/M + F > M\n"
            "   X:2699521-154931043     M/M + F > M/F\n"
            "   Y:1-59373566            .   + F > F\n"
            "   MT:1-16569              M   + F > M\n"
            "\n"
            "   chrX:1-60000            M/M + F > M\n"
            "   chrX:1-60000            M/M + F > M/F\n"
            "   chrX:2699521-154931043  M/M + F > M\n"
            "   chrX:2699521-154931043  M/M + F > M/F\n"
            "   chrY:1-59373566         .   + F > F\n"
            "   chrM:1-16569            M   + F > M\n"
    },
    { .alias = "GRCh38",
      .about = "Human Genome reference assembly GRCh38 / hg38, both chr naming conventions",
      .rules =
            "   X:1-9999                M/M + F > M\n"
            "   X:1-9999                M/M + F > M/F\n"
            "   X:2781480-155701381     M/M + F > M\n"
            "   X:2781480-155701381     M/M + F > M/F\n"
            "   Y:1-57227415            .   + F > F\n"
            "   MT:1-16569              M   + F > M\n"
            "\n"
            "   chrX:1-9999             M/M + F > M\n"
            "   chrX:1-9999             M/M + F > M/F\n"
            "   chrX:2781480-155701381  M/M + F > M\n"
            "   chrX:2781480-155701381  M/M + F > M/F\n"
            "   chrY:1-57227415         .   + F > F\n"
            "   chrM:1-16569            M   + F > M\n"
    },
    {
        .alias = NULL,
        .about = NULL,
        .rules = NULL,
    }
};


const char *usage(void)
{
    return
        "\n"
        "About: Count Mendelian consistent / inconsistent genotypes. Note that this plugin is DEPRECATED and\n"
        "       will not be supported in the future. Please use the newer plugin +mendelian2 instead.\n"
        "Usage: bcftools +mendelian [Options]\n"
        "Options:\n"
        "   -c, --count                     Count the number of consistent sites [DEPRECATED, use `-m c` instead]\n"
        "   -d, --delete                    Delete inconsistent genotypes (set to \"./.\") [DEPRECATED, use `-m d` instead]\n"
        "   -l, --list [+x]                 List consistent (+) or inconsistent (x) sites [DEPRECATED, use `-m +` or `-m x` instead]\n"
        "   -m, --mode [+acdux]             Output mode (the default is `-m c`):\n"
        "                                       + .. list consistent sites\n"
        "                                       a .. add INFO/MERR annotation with the number of inconsistent samples\n"
        "                                       c .. print counts, a text summary with the number of errors per trio\n"
        "                                       d .. delete inconsistent genotypes (set to \"./.\")\n"
        "                                       u .. list uninformative sites\n"
        "                                       x .. list inconsistent sites\n"
        "   -o, --output FILE               Write output to a file [standard output]\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "   -r, --rules ASSEMBLY[?]         Predefined rules, 'list' to print available settings, append '?' for details\n"
        "   -R, --rules-file FILE           Inheritance rules, see example below\n"
        "   -t, --trio M,F,C                Names of mother, father and the child\n"
        "   -T, --trio-file FILE            List of trios, one per line (mother,father,child)\n"
        "   -p, --ped FILE                  PED file\n"
        "\n"
        "Example:\n"
        "   # Default inheritance patterns, override with -r\n"
        "   #   region  maternal_ploidy + paternal > offspring\n"
        "   X:1-60000            M/M + F > M\n"
        "   X:1-60000            M/M + F > M/F\n"
        "   X:2699521-154931043  M/M + F > M\n"
        "   X:2699521-154931043  M/M + F > M/F\n"
        "   Y:1-59373566         .   + F > F\n"
        "   MT:1-16569           M   + F > M\n"
        "\n"
        "   bcftools +mendelian in.vcf -t Mother,Father,Child -c\n"
        "\n";
}

regidx_t *init_rules(args_t *args, char *alias)
{
    const rules_predef_t *rules = rules_predefs;
    if ( !alias ) alias = "GRCh37";

    int detailed = 0, len = strlen(alias);
    if ( alias[len-1]=='?' ) { detailed = 1; alias[len-1] = 0; }

    while ( rules->alias && strcasecmp(alias,rules->alias) ) rules++;

    if ( !rules->alias )
    {
        fprintf(stderr,"\nPRE-DEFINED INHERITANCE RULES\n\n");
        fprintf(stderr," * Columns are: CHROM:BEG-END MATERNAL_PLOIDY + PATERNAL_PLOIDY > OFFSPRING\n");
        fprintf(stderr," * Coordinates are 1-based inclusive.\n\n");
        rules = rules_predefs;
        while ( rules->alias )
        {
            fprintf(stderr,"%s\n   .. %s\n\n", rules->alias,rules->about);
            if ( detailed )
                fprintf(stderr,"%s\n", rules->rules);
            rules++;
        }
        fprintf(stderr,"Run as --rules <alias> (e.g. --rules GRCh37).\n");
        fprintf(stderr,"To see the detailed ploidy definition, append a question mark (e.g. --rules GRCh37?).\n");
        fprintf(stderr,"\n");
        exit(-1);
    }
    else if ( detailed )
    {
        fprintf(stderr,"%s", rules->rules);
        exit(-1);
    }
    return regidx_init_string(rules->rules, parse_rules, NULL, sizeof(rule_t), &args);
}

static int parse_rules(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    // e.g. "Y:1-59373566        .   + F > . # daugther"

    // eat any leading spaces
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) return -1;      // skip empty lines

    // chromosome name, beg, end
    char *tmp, *se = ss;
    while ( se[1] && !isspace(se[1]) ) se++;
    while ( se > ss && isdigit(*se) ) se--;
    if ( *se!='-' ) error("Could not parse the region: %s\n", line);
    *end = strtol(se+1, &tmp, 10) - 1;
    if ( tmp==se+1 ) error("Could not parse the region:%s\n",line);
    while ( se > ss && *se!=':' ) se--;
    *beg = strtol(se+1, &tmp, 10) - 1;
    if ( tmp==se+1 ) error("Could not parse the region:%s\n",line);

    *chr_beg = ss;
    *chr_end = se-1;

    // skip region
    while ( *ss && !isspace(*ss) ) ss++;
    while ( *ss && isspace(*ss) ) ss++;

    rule_t *rule = (rule_t*) payload;
    memset(rule, 0, sizeof(rule_t));

    // maternal ploidy
    se = ss;
    while ( *se && !isspace(*se) ) se++;
    int err = 0;
    if ( se - ss == 1 )
    {
        if ( *ss=='M' ) rule->mpl = 1;
        else if ( *ss=='.' ) rule->mpl = 0;
        else err = 1;
    }
    else if ( se - ss == 3 )
    {
        if ( !strncmp(ss,"M/M",3) ) rule->mpl = 2;
        else err = 1;
    }
    else err = 1;
    if ( err ) error("Could not parse the maternal ploidy, only \"M\", \"M/M\" and \".\" currently supported: %s\n",line);

    // skip "+"
    while ( *se && isspace(*se) ) se++;
    if ( *se != '+' ) error("Could not parse the line: %s\n",line);
    se++;
    while ( *se && isspace(*se) ) se++;

    // paternal ploidy
    ss = se;
    while ( *se && !isspace(*se) ) se++;
    if ( se - ss == 1 )
    {
        if ( *ss=='F' ) rule->fpl = 1;
        else err = 1;
    }
    else err = 1;
    if ( err ) error("Could not parse the paternal ploidy, only \"F\" is currently supported: %s [%s]\n",line, ss);

    // skip ">"
    while ( *se && isspace(*se) ) se++;
    if ( *se != '>' ) error("Could not parse the line: %s\n",line);
    se++;
    while ( *se && isspace(*se) ) se++;

    // ploidy in offspring
    ss = se;
    while ( *se && !isspace(*se) ) se++;
    if ( se - ss == 3 )
    {
        if ( !strncmp(ss,"M/F",3) ) { rule->cpl = 2; rule->fal = 1; rule->mal = 1; }
        else err = 1;
    }
    else if ( se - ss == 1 )
    {
        if ( *ss=='F' ) { rule->cpl = 1; rule->fal = 1; }
        else if ( *ss=='M' ) { rule->cpl = 1; rule->mal = 1; }
        else err = 1;
    }
    else err = 1;
    if ( err ) error("Could not parse the offspring's ploidy, only \"M\", \"F\" or \"M/F\" is currently supported: %s\n",line);

    return 0;
}

void parse_ped(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL;
    do
    {
        int ncols = ksplit_core(str.s,0,&moff,&off);
        if ( ncols<4 ) error("Could not parse the ped file: %s\n", str.s);

        int ifather = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[2]]);
        int imother = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[3]]);
        int ichild = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[1]]);

        // The code in process() makes an attempt to work with partial families,
        // the support is not complete though and can lead to core dumps. Therefore
        // enforcing full trios for now.
        // if ( ( ifather<0 && imother<0 ) || ichild<0 ) continue;
        if ( ifather<0 || imother<0 || ichild<0 ) continue;

        args->ntrios++;
        hts_expand0(trio_t,args->ntrios,args->mtrios,args->trios);
        trio_t *trios = &args->trios[args->ntrios-1];
        trios->ifather = ifather;
        trios->imother = imother;
        trios->ichild  = ichild;

    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );
    if ( !args->ntrios ) error("No complete trios found in the PED and VCF\n");

    free(str.s);
    free(off);
    hts_close(fp);
}

int run(int argc, char **argv)
{
    char *trio_samples = NULL, *trio_file = NULL, *ped_fname = NULL, *rules_fname = NULL, *rules_string = NULL;
    memset(&args,0,sizeof(args_t));
    args.mode = 0;
    args.output_fname = "-";
    args.clevel = -1;

    static struct option loptions[] =
    {
        {"trio",1,0,'t'},
        {"trio-file",1,0,'T'},
        {"ped",1,0,'p'},
        {"delete",0,0,'d'},
        {"list",1,0,'l'},
        {"mode",1,0,'m'},
        {"count",0,0,'c'},
        {"rules",1,0,'r'},
        {"rules-file",1,0,'R'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {0,0,0,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "?ht:T:p:l:m:cdr:R:o:O:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'o': args.output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args.output_type = FT_BCF_GZ; break;
                          case 'u': args.output_type = FT_BCF; break;
                          case 'z': args.output_type = FT_VCF_GZ; break;
                          case 'v': args.output_type = FT_VCF; break;
                          default:
                          {
                              args.clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args.clevel<0 || args.clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      };
                      if ( optarg[1] )
                      {
                          args.clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args.clevel<0 || args.clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case 'R': rules_fname = optarg; break;
            case 'r': rules_string = optarg; break;
            case 'd':
                args.mode |= MODE_DELETE;
                fprintf(stderr,"Warning: -d will be deprecated, please use `-m d` instead.\n");
                break;
            case 'c':
                args.mode |= MODE_COUNT;
                fprintf(stderr,"Warning: -c will be deprecated, please use `-m c` instead.\n");
                break;
            case 'l':
                if ( !strcmp("+",optarg) ) args.mode |= MODE_LIST_GOOD;
                else if ( !strcmp("x",optarg) ) args.mode |= MODE_LIST_BAD;
                else error("The argument not recognised: --list %s\n", optarg);
                fprintf(stderr,"Warning: -l will be deprecated, please use -m instead.\n");
                break;
            case 'm':
                if ( !strcmp("+",optarg) ) args.mode |= MODE_LIST_GOOD;
                else if ( !strcmp("x",optarg) ) args.mode |= MODE_LIST_BAD;
                else if ( !strcmp("a",optarg) ) args.mode |= MODE_ANNOTATE;
                else if ( !strcmp("d",optarg) ) args.mode |= MODE_DELETE;
                else if ( !strcmp("c",optarg) ) args.mode |= MODE_COUNT;
                else if ( !strcmp("u",optarg) ) args.mode |= MODE_LIST_SKIP;
                else error("The argument not recognised: --mode %s\n", optarg);
                break;
            case 't': trio_samples = optarg; break;
            case 'T': trio_file = optarg; break;
            case 'p': ped_fname = optarg; break;
            case 'h':
            case '?':
            default: error("%s",usage()); break;
        }
    }
    if ( rules_fname )
        args.rules = regidx_init(rules_fname, parse_rules, NULL, sizeof(rule_t), &args);
    else
        args.rules = init_rules(&args, rules_string);
    if ( !args.rules ) return -1;
    args.itr     = regitr_init(args.rules);
    args.itr_ori = regitr_init(args.rules);

    char *fname = NULL;
    if ( optind>=argc || argv[optind][0]=='-' )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else error("%s",usage());
    }
    else
        fname = argv[optind];

    if ( !trio_samples && !trio_file && !ped_fname ) error("Expected the -t/T or -p option\n");
    if ( !args.mode ) args.mode = MODE_COUNT;
    if ( args.mode&MODE_DELETE && !(args.mode&(MODE_LIST_GOOD|MODE_LIST_BAD|MODE_LIST_SKIP)) ) args.mode |= MODE_LIST_GOOD|MODE_LIST_BAD|MODE_LIST_SKIP;
    if ( args.mode&MODE_ANNOTATE && !(args.mode&(MODE_LIST_GOOD|MODE_LIST_BAD|MODE_LIST_SKIP)) ) args.mode |= MODE_LIST_GOOD|MODE_LIST_BAD|MODE_LIST_SKIP;

    FILE *log_fh = stderr;
    if ( args.mode==MODE_COUNT )
    {
        log_fh = strcmp("-",args.output_fname) ? fopen(args.output_fname,"w") : stdout;
        if ( !log_fh ) error("Error: cannot write to %s\n", args.output_fname);
    }

    args.sr = bcf_sr_init();
    if ( !bcf_sr_add_reader(args.sr, fname) ) error("Failed to read from %s: %s\n", !strcmp("-",fname)?"standard input":fname,bcf_sr_strerror(args.sr->errnum));
    args.hdr = bcf_sr_get_header(args.sr, 0);
    if ( args.mode!=MODE_COUNT )
    {
        char wmode[8];
        set_wmode(wmode,args.output_type,args.output_fname,args.clevel);
        args.out_fh = hts_open(args.output_fname ? args.output_fname : "-", wmode);
        if ( args.out_fh == NULL ) error("Can't write to \"%s\": %s\n", args.output_fname, strerror(errno));
        if ( args.mode&MODE_ANNOTATE )
            bcf_hdr_append(args.hdr, "##INFO=<ID=MERR,Number=1,Type=Integer,Description=\"Number of trios with Mendelian errors\">");
        if ( bcf_hdr_write(args.out_fh, args.hdr)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args.output_fname);
    }

    int i, n = 0;
    char **list;
    if ( trio_samples )
    {
        args.ntrios = 1;
        args.trios = (trio_t*) calloc(1,sizeof(trio_t));
        list = hts_readlist(trio_samples, 0, &n);
        if ( n!=3 ) error("Expected three sample names with -t\n");
        args.trios[0].imother = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[0]);
        args.trios[0].ifather = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[1]);
        args.trios[0].ichild  = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[2]);
        if ( args.trios[0].imother<0 ) error("The sample is not present in the VCF: %s\n",list[0]);
        if ( args.trios[0].ifather<0 ) error("The sample is not present in the VCF: %s\n",list[1]);
        if ( args.trios[0].ichild<0 )  error("The sample is not present in the VCF: %s\n",list[2]);
        for (i=0; i<n; i++) free(list[i]);
        free(list);
    }
    if ( trio_file )
    {
        list = hts_readlist(trio_file, 1, &n);
        if ( !list ) error("Error: could not read file %s\n",trio_file);
        args.ntrios = n;
        args.trios = (trio_t*) calloc(n,sizeof(trio_t));
        for (i=0; i<n; i++)
        {
            char *ss = list[i], *se;
            se = strchr(ss, ',');
            if ( !se ) error("Could not parse %s: %s\n",trio_file, ss);
            *se = 0;
            args.trios[i].imother = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, ss);
            if ( args.trios[i].imother<0 ) error("No such sample: \"%s\"\n", ss);
            ss = ++se;
            se = strchr(ss, ',');
            if ( !se ) error("Could not parse %s\n",trio_file);
            *se = 0;
            args.trios[i].ifather = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, ss);
            if ( args.trios[i].ifather<0 ) error("No such sample: \"%s\"\n", ss);
            ss = ++se;
            if ( *ss=='\0' ) error("Could not parse %s\n",trio_file);
            args.trios[i].ichild = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, ss);
            if ( args.trios[i].ichild<0 ) error("No such sample: \"%s\"\n", ss);
            free(list[i]);
        }
        free(list);
    }
    if ( ped_fname ) parse_ped(&args, ped_fname);

    while ( bcf_sr_next_line(args.sr) )
    {
        bcf1_t *line = bcf_sr_get_line(args.sr,0);
        line = process(line);
        if ( line )
        {
            if ( line->errcode ) error("TODO: Unchecked error (%d), exiting\n",line->errcode);
            if ( args.out_fh && bcf_write1(args.out_fh, args.hdr, line)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args.output_fname);
        }
    }
    if ( args.out_fh && hts_close(args.out_fh)!=0 ) error("Error: close failed\n");

    if ( args.mode & MODE_COUNT )
    {
        fprintf(log_fh,"# [1]nOK\t[2]nBad\t[3]nSkipped\t[4]Trio (mother,father,child)\n");
        for (i=0; i<args.ntrios; i++)
        {
            trio_t *trio = &args.trios[i];
            fprintf(log_fh,"%d\t%d\t%d\t%s,%s,%s\n",
                    trio->nok,trio->nbad,args.nrec-(trio->nok+trio->nbad),
                    bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, trio->imother),
                    bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, trio->ifather),
                    bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, trio->ichild)
                   );
        }
    }
    if ( log_fh!=stderr && log_fh!=stdout && fclose(log_fh) ) error("Error: close failed for %s\n", args.output_fname);

    free(args.gt_arr);
    free(args.trios);
    regitr_destroy(args.itr);
    regitr_destroy(args.itr_ori);
    regidx_destroy(args.rules);
    bcf_sr_destroy(args.sr);
    return 0;
}

static void warn_ploidy(bcf1_t *rec)
{
    static int warned = 0;
    if ( warned ) return;
    fprintf(stderr,"Incorrect ploidy at %s:%"PRId64", skipping the trio. (This warning is printed only once.)\n", bcf_seqname(args.hdr,rec),(int64_t) rec->pos+1);
    warned = 1;
}

bcf1_t *process(bcf1_t *rec)
{
    bcf1_t *dflt = args.mode&MODE_LIST_SKIP ? rec : NULL;
    args.nrec++;

    if ( rec->n_allele > 63 ) return dflt;      // we use 64bit bitmask below

    int ngt = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.ngt_arr);
    if ( ngt<0 ) return dflt;
    if ( ngt!=2*bcf_hdr_nsamples(args.hdr) && ngt!=bcf_hdr_nsamples(args.hdr) ) return dflt;
    ngt /= bcf_hdr_nsamples(args.hdr);

    int itr_set = regidx_overlap(args.rules, bcf_seqname(args.hdr,rec),rec->pos,rec->pos, args.itr_ori);

    int i, nbad = 0, ngood = 0, needs_update = 0;
    for (i=0; i<args.ntrios; i++)
    {
        int32_t a,b,c,d,e,f;
        trio_t *trio = &args.trios[i];

        if ( trio->imother<0 )
        {
            a = bcf_gt_missing;
            b = bcf_int32_vector_end;
        }
        else
        {
            a = args.gt_arr[ngt*trio->imother];
            b = ngt==2 ? args.gt_arr[ngt*trio->imother+1] : bcf_int32_vector_end;
        }
        if ( trio->ifather<0 )
        {
            c = bcf_gt_missing;
            d = bcf_int32_vector_end;
        }
        else
        {
          c = args.gt_arr[ngt*trio->ifather];
          d = ngt==2 ? args.gt_arr[ngt*trio->ifather+1] : bcf_int32_vector_end;
        }
        e = args.gt_arr[ngt*trio->ichild];
        f = ngt==2 ? args.gt_arr[ngt*trio->ichild+1] : bcf_int32_vector_end;

        // skip sites with missing data in child
        if ( bcf_gt_is_missing(e) || bcf_gt_is_missing(f) ) continue;

        uint64_t mother = 0, father = 0,child1,child2;

        int is_ok = 0;
        if ( !itr_set )
        {
            if ( f==bcf_int32_vector_end ) { warn_ploidy(rec); continue; }

            // All M,F,C genotypes are diploid. Missing data are considered consistent.
            child1 = 1<<bcf_gt_allele(e);
            child2 = 1<<bcf_gt_allele(f);
            mother  = bcf_gt_is_missing(a) ? child1|child2 : 1<<bcf_gt_allele(a);
            mother |= bcf_gt_is_missing(b) || b==bcf_int32_vector_end ? child1|child2 : 1<<bcf_gt_allele(b);
            father  = bcf_gt_is_missing(c) ? child1|child2 : 1<<bcf_gt_allele(c);
            father |= bcf_gt_is_missing(d) || d==bcf_int32_vector_end ? child1|child2 : 1<<bcf_gt_allele(d);

            if ( (mother&child1 && father&child2) || (mother&child2 && father&child1) ) is_ok = 1;
        }
        else
        {
            child1  = 1<<bcf_gt_allele(e);
            child2  = bcf_gt_is_missing(f) || f==bcf_int32_vector_end ? 0 : 1<<bcf_gt_allele(f);
            mother |= bcf_gt_is_missing(a) ? 0 : 1<<bcf_gt_allele(a);
            mother |= bcf_gt_is_missing(b) || b==bcf_int32_vector_end ? 0 : 1<<bcf_gt_allele(b);
            father |= bcf_gt_is_missing(c) ? 0 : 1<<bcf_gt_allele(c);
            father |= bcf_gt_is_missing(d) || d==bcf_int32_vector_end ? 0 : 1<<bcf_gt_allele(d);

            regitr_copy(args.itr, args.itr_ori);
            while ( !is_ok && regitr_overlap(args.itr) )
            {
                rule_t *rule = &regitr_payload(args.itr,rule_t);
                if ( child1 && child2 )
                {
                    if ( !rule->mal || !rule->fal ) continue;   // wrong rule (haploid), but this is a diploid GT
                    if ( !mother ) mother = child1|child2;
                    if ( !father ) father = child1|child2;
                    if ( (mother&child1 && father&child2) || (mother&child2 && father&child1) ) is_ok = 1;
                    continue;
                }
                if ( rule->mal )
                {
                    if ( mother && !(child1&mother) ) continue;
                }
                if ( rule->fal )
                {
                    if ( father && !(child1&father) ) continue;
                }
                is_ok = 1;
            }
        }
        if ( is_ok )
        {
            trio->nok++;
            ngood++;
        }
        else
        {
            trio->nbad++;
            nbad++;
            if ( args.mode&MODE_DELETE )
            {
                args.gt_arr[ngt*trio->imother] = bcf_gt_missing;
                if ( b!=bcf_int32_vector_end ) args.gt_arr[ngt*trio->imother+1] = bcf_gt_missing; // should be always true
                args.gt_arr[ngt*trio->ifather] = bcf_gt_missing;
                if ( d!=bcf_int32_vector_end ) args.gt_arr[ngt*trio->ifather+1] = bcf_gt_missing;
                args.gt_arr[ngt*trio->ichild] = bcf_gt_missing;
                if ( f!=bcf_int32_vector_end ) args.gt_arr[ngt*trio->ichild+1]  = bcf_gt_missing;
                needs_update = 1;
            }
        }
    }

    if ( needs_update && bcf_update_genotypes(args.hdr,rec,args.gt_arr,ngt*bcf_hdr_nsamples(args.hdr)) )
        error("Could not update GT field at %s:%"PRId64"\n", bcf_seqname(args.hdr,rec),(int64_t) rec->pos+1);

    if ( args.mode&MODE_ANNOTATE ) bcf_update_info_int32(args.hdr, rec, "MERR", &nbad, 1);
    if ( args.mode&MODE_LIST_GOOD && ngood ) return rec;
    if ( args.mode&MODE_LIST_BAD && nbad ) return rec;
    if ( args.mode&MODE_LIST_SKIP && !ngood && !nbad ) return rec;

    return NULL;
}
