/* The MIT License

   Copyright (c) 2025 Gekkonid Scientific Pty. Ltd.

   Author: Kevin Murray <kevin@gekkonid.com>

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

#include <getopt.h>
#include <htslib/kbitset.h>
#include <htslib/vcf.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

typedef struct {
    float minAC;
    float maxAC;
    float minAF;
    float maxAF;
    float minAD;
    float maxAD;
    uint64_t n_sites_trimmed;
    uint64_t n_sites;
    int has_ac_filter;
    int has_af_filter;
    int has_ad_filter;
    bcf_hdr_t* hdr_in;
    bcf_hdr_t* hdr_out;
    kbitset_t* warn_tag;
} pa_cfg_t;

static pa_cfg_t cfg = {
    .minAC = -1,
    .maxAC = -1,
    .minAF = -1,
    .maxAF = -1,
    .minAD = -1,
    .maxAD = -1,
    .n_sites_trimmed = 0,
    .n_sites = 0,
    .has_ac_filter = 0,
    .has_af_filter = 0,
    .has_ad_filter = 0,
    .hdr_in = NULL,
    .hdr_out = NULL
};

const char* about(void)
{
    return "Prune multialleic sites by AC, AF, or AD, replacing pruned allels with missing genotypes.\n"
           "Useful in cases where the Nth allele of a multialleic is exceptionally rare or likely erroeous.\n";
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About: %s", about());
    fprintf(stderr, "Usage: bcftools +prune-alleles [options] <in.vcf.gz>\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "    -a, --min-ac <float>       minimum allele count (inclusive)\n");
    fprintf(stderr, "    -A, --max-ac <float>       maximum allele count (inclusive)\n");
    fprintf(stderr, "    -f, --min-af <float>       minimum allele frequency (inclusive)\n");
    fprintf(stderr, "    -F, --max-af <float>       maximum allele frequency (inclusive)\n");
    fprintf(stderr, "    -d, --min-ad <float>       minimum allele depth (inclusive)\n");
    fprintf(stderr, "    -D, --max-ad <float>       maximum allele depth (inclusive)\n");
    fprintf(stderr, "\n");
    exit(1);
}

/*
    Called once at startup, it initializes local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(int argc, char** argv, bcf_hdr_t* in, bcf_hdr_t* out)
{
    int c;
    static struct option loptions[] = {
        { "min-ac", required_argument, 0, 'a' },
        { "max-ac", required_argument, 0, 'A' },
        { "min-af", required_argument, 0, 'f' },
        { "max-af", required_argument, 0, 'F' },
        { "min-ad", required_argument, 0, 'd' },
        { "max-ad", required_argument, 0, 'D' },
        { 0, 0, 0, 0 }
    };

    while ((c = getopt_long(argc, argv, "a:A:f:F:d:D:h?", loptions, NULL)) >= 0) {
        switch (c) {
        case 'a':
            cfg.minAC = atof(optarg);
            cfg.has_ac_filter = 1;
            break;
        case 'A':
            cfg.maxAC = atof(optarg);
            cfg.has_ac_filter = 1;
            break;
        case 'f':
            cfg.minAF = atof(optarg);
            cfg.has_af_filter = 1;
            break;
        case 'F':
            cfg.maxAF = atof(optarg);
            cfg.has_af_filter = 1;
            break;
        case 'd':
            cfg.minAD = atof(optarg);
            cfg.has_ad_filter = 1;
            break;
        case 'D':
            cfg.maxAD = atof(optarg);
            cfg.has_ad_filter = 1;
            break;
        case 'h':
        case '?':
        default:
            usage();
            break;
        }
    }

    // Check if at least one filter is set
    if (!cfg.has_ac_filter && !cfg.has_af_filter && !cfg.has_ad_filter) {
        fprintf(stderr, "Error: At least one filter must be specified\n");
        usage();
    }

    // Ensure min/max values are consistent
    if (cfg.has_ac_filter && cfg.minAC > 0 && cfg.maxAC > 0 && cfg.minAC > cfg.maxAC) {
        fprintf(stderr, "Error: --min-ac cannot be greater than --max-ac\n");
        exit(1);
    }
    if (cfg.has_af_filter && cfg.minAF > 0 && cfg.maxAF > 0 && cfg.minAF > cfg.maxAF) {
        fprintf(stderr, "Error: --min-af cannot be greater than --max-af\n");
        exit(1);
    }
    if (cfg.has_ad_filter && cfg.minAD > 0 && cfg.maxAD > 0 && cfg.minAD > cfg.maxAD) {
        fprintf(stderr, "Error: --min-ad cannot be greater than --max-ad\n");
        exit(1);
    }

    // Ensure we have required header
    cfg.hdr_in = in;
    cfg.hdr_out = out;
    if (!cfg.hdr_in) {
        fprintf(stderr, "Error: Can't read V/BCF header. Are you providing data?\n");
        usage();
        exit(1);
    }
    int res;
    if (cfg.has_ac_filter) {
        res = bcf_hdr_id2int(cfg.hdr_in, BCF_DT_ID, "AC");
        if (res < 0 || !bcf_hdr_idinfo_exists(cfg.hdr_in, BCF_HL_INFO, res)) {
            fprintf(stderr, "Error: --min-ac/--max-ac require pre-filled INFO/AC, use bcftools +fill-tags :)\n");
            exit(1);
        }
    }
    if (cfg.has_af_filter) {
        res = bcf_hdr_id2int(cfg.hdr_in, BCF_DT_ID, "AF");
        if (res < 0 || !bcf_hdr_idinfo_exists(cfg.hdr_in, BCF_HL_INFO, res)) {
            fprintf(stderr, "Error: --min-af/--max-af require pre-filled INFO/AF, use bcftools +fill-tags :)\n");
            exit(1);
        }
    }
    if (cfg.has_ad_filter) {
        res = bcf_hdr_id2int(cfg.hdr_in, BCF_DT_ID, "AD");
        if (res < 0 || !bcf_hdr_idinfo_exists(cfg.hdr_in, BCF_HL_INFO, res)) {
            fprintf(stderr, "Error: --min-ad/--max-ad require pre-filled INFO/AD, use bcftools +fill-tags :)\n");
            exit(1);
        }
    }

    cfg.warn_tag = kbs_init((cfg.hdr_in)->n[BCF_DT_ID]);
    return 0;
}

/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t* process(bcf1_t* rec)
{
    cfg.n_sites++;

    // Skip if not multiallelic
    if (rec->n_allele <= 2)
        return rec;

    bcf_unpack(rec, BCF_UN_ALL);

    // Initialize arrays to track which alleles to keep
    kbitset_t* keep_allele = kbs_init(rec->n_allele);
    kbs_insert(keep_allele, 0); // Always keep the reference allele

    int i, j, k, ret, n_ac = 0, n_af = 0, n_ad = 0;
    int32_t *ac = NULL, *ad = NULL;
    float* af = NULL;

    // Get AC values if needed
    if (cfg.has_ac_filter) {
        ret = bcf_get_info_int32(cfg.hdr_in, rec, "AC", &ac, &n_ac);
        if (ret <= 0 || n_ac != rec->n_allele - 1) {
            if (ac)
                free(ac);
            ac = NULL;
            n_ac = 0;
        }
    }

    // Get AF values if needed
    if (cfg.has_af_filter) {
        ret = bcf_get_info_float(cfg.hdr_in, rec, "AF", &af, &n_af);
        if (ret <= 0 || n_af != rec->n_allele - 1) {
            if (af)
                free(af);
            af = NULL;
            n_af = 0;
        }
    }

    // Get AD values if needed
    if (cfg.has_ad_filter) {
        ret = bcf_get_info_int32(cfg.hdr_in, rec, "AD", &ad, &n_ad);
        if (ret <= 0 || n_ad != rec->n_allele) {
            if (ad)
                free(ad);
            ad = NULL;
            n_ad = 0;
        }
    }

    // Check if we have any data to filter on
    if ((!ac || n_ac == 0) && (!af || n_af == 0) && (!ad || n_ad == 0)) {
        kbs_destroy(keep_allele);
        if (ac)
            free(ac);
        if (af)
            free(af);
        if (ad)
            free(ad);
        return rec;
    }

    // Determine which alleles to keep
    int any_filtered = 0;
    for (i = 1; i < rec->n_allele; i++) {
        int keep = 1;

        // Check AC
        if (ac && n_ac > 0 && i - 1 < n_ac) {
            if ((cfg.minAC > 0 && ac[i - 1] < cfg.minAC) || (cfg.maxAC > 0 && ac[i - 1] > cfg.maxAC)) {
                keep = 0;
            }
        }

        // Check AF
        if (keep && af && n_af > 0 && i - 1 < n_af) {
            if ((cfg.minAF > 0 && af[i - 1] < cfg.minAF) || (cfg.maxAF > 0 && af[i - 1] > cfg.maxAF)) {
                keep = 0;
            }
        }

        // Check AD
        if (keep && ad && n_ad > 0 && i < n_ad) {
            if ((cfg.minAD > 0 && ad[i] < cfg.minAD) || (cfg.maxAD > 0 && ad[i] > cfg.maxAD)) {
                keep = 0;
            }
        }

        if (keep) {
            kbs_insert(keep_allele, i);
        } else {
            any_filtered = 1;
        }
    }

    // If no alleles were filtered, return the record as is
    if (!any_filtered) {
        kbs_destroy(keep_allele);
        if (ac)
            free(ac);
        if (af)
            free(af);
        if (ad)
            free(ad);
        return rec;
    }

    // Count how many alleles we're keeping
    int n_keep = 0;
    int n_allele_orig = rec->n_allele;
    for (i = 0; i < rec->n_allele; i++) {
        if (kbs_exists(keep_allele, i))
            n_keep++;
    }

    if (n_keep <= 1) {
        kbs_destroy(keep_allele);
        if (ac)
            free(ac);
        if (af)
            free(af);
        if (ad)
            free(ad);
        return NULL;
    }

    // Create new arrays for alleles and INFO fields
    char** alleles = (char**)malloc(n_keep * sizeof(char*));
    int n_new = 0;

    // Map old allele indices to new ones
    int* old2new = (int*)malloc(rec->n_allele * sizeof(int));
    for (i = 0; i < rec->n_allele; i++) {
        if (kbs_exists(keep_allele, i)) {
            alleles[n_new] = rec->d.allele[i];
            old2new[i] = n_new++;
        } else {
            old2new[i] = -1;
        }
    }

    // Update all INFO fields with Number=A (one value per alternate allele) or Number=R (one value per allele)
    for (j = 0; j < rec->n_info; j++) {
        bcf_info_t* info = &rec->d.info[j];
        const char* key = cfg.hdr_in->id[BCF_DT_ID][info->key].key;
        int hdr_id = info->key;

        // Use bcf_hdr_id2length to determine if this is a Number=A or Number=R field
        int id_len = bcf_hdr_id2length(cfg.hdr_in, BCF_HL_INFO, hdr_id);
        if (!(id_len == BCF_VL_A || id_len == BCF_VL_R))
            continue;

        // Count how many alleles we're keeping (total for R, alt-only for A)
        int n_values_to_keep = (id_len == BCF_VL_R) ? n_keep : (n_keep - 1);

        // Handle different types of INFO fields
        if (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32) {
            int n_val = 0;
            int32_t* vals = NULL;
            int ret = bcf_get_info_int32(cfg.hdr_in, rec, key, &vals, &n_val);
            int expected_values = (id_len == BCF_VL_R) ? n_allele_orig : (n_allele_orig - 1);

            if (ret > 0 && n_val == expected_values) {
                int32_t* new_vals = (int32_t*)malloc(n_values_to_keep * sizeof(int32_t));
                int new_idx = 0;
                int start_idx = (id_len == BCF_VL_R) ? 0 : 1;

                for (i = start_idx; i < n_allele_orig; i++) {
                    if (kbs_exists(keep_allele, i)) {
                        new_vals[new_idx++] = vals[i - start_idx];
                    }
                }

                bcf_update_info_int32(cfg.hdr_in, rec, key, new_vals, n_values_to_keep);
                free(new_vals);
            } else if (!kbs_exists(cfg.warn_tag, hdr_id)) {
                kbs_insert(cfg.warn_tag, hdr_id);
                fprintf(stderr, "WARNING: got unexpected number of value from INFO/%s, so leaving it intact: (got %d, expected %d)\n", key, n_val, expected_values);
            }
            if (vals)
                free(vals);
        } else if (info->type == BCF_BT_FLOAT) {
            int n_val = 0;
            float* vals = NULL;
            int ret = bcf_get_info_float(cfg.hdr_in, rec, key, &vals, &n_val);
            int expected_values = (id_len == BCF_VL_R) ? n_allele_orig : (n_allele_orig - 1);

            if (ret > 0 && n_val == expected_values) {
                float* new_vals = (float*)malloc(n_values_to_keep * sizeof(float));
                int new_idx = 0;
                int start_idx = (id_len == BCF_VL_R) ? 0 : 1;

                for (i = start_idx; i < n_allele_orig; i++) {
                    if (kbs_exists(keep_allele, i)) {
                        new_vals[new_idx++] = vals[i - start_idx];
                    }
                }

                bcf_update_info_float(cfg.hdr_in, rec, key, new_vals, n_values_to_keep);
                free(new_vals);
            } else if (!kbs_exists(cfg.warn_tag, hdr_id)) {
                kbs_insert(cfg.warn_tag, hdr_id);
                fprintf(stderr, "WARNING: got unexpected number of value from INFO/%s, so leaving it intact: (got %d, expected %d)\n", key, n_val, expected_values);
            }
            if (vals)
                free(vals);
        } else if (info->type == BCF_BT_CHAR) {
            int n_val = 0;
            char* vals = NULL;
            int ret = bcf_get_info_string(cfg.hdr_in, rec, key, &vals, &n_val);
            if (ret > 0) {
                // String fields are more complex, we'd need to parse the comma-separated values
                // For now, we'll just remove the field as it's less common
                bcf_update_info_string(cfg.hdr_in, rec, key, NULL);
                if (!kbs_exists(cfg.warn_tag, hdr_id)) {
                    kbs_insert(cfg.warn_tag, hdr_id);
                    fprintf(stderr, "WARNING: Removing string multiallelic field INFO/%s as it will be inaccurate\n", key);
                }
            }
            if (vals)
                free(vals);
        }
    }

    // Update all FORMAT fields with Number=R (one value per allele) or Number=A (one value per alt allele)
    int n_samples = bcf_hdr_nsamples(cfg.hdr_in);

    // Get all FORMAT fields
    for (j = 0; j < rec->n_fmt; j++) {
        bcf_fmt_t* fmt = &rec->d.fmt[j];
        const char* key = cfg.hdr_in->id[BCF_DT_ID][fmt->id].key;
        int hdr_id = fmt->id;

        // Use bcf_hdr_id2length to determine if this is a Number=A or Number=R field
        int id_len = bcf_hdr_id2length(cfg.hdr_in, BCF_HL_FMT, hdr_id);
        if (!(id_len == BCF_VL_A || id_len == BCF_VL_R))
            continue;
        int field_type = bcf_hdr_id2type(cfg.hdr_in, BCF_HL_FMT, hdr_id);

        if (field_type == BCF_HT_INT) {
            int n_val = 0;
            int32_t* vals = NULL;
            ret = bcf_get_format_int32(cfg.hdr_in, rec, key, &vals, &n_val);
            if (ret > 0) {
                int n_values_per_sample = ret / n_samples;
                int expected_values = (id_len == BCF_VL_R) ? n_allele_orig : (n_allele_orig - 1);

                if (n_values_per_sample == expected_values) {
                    int new_n_values = (id_len == BCF_VL_R) ? n_keep : (n_keep - 1);
                    int32_t* new_vals = (int32_t*)malloc(n_samples * new_n_values * sizeof(int32_t));

                    for (i = 0; i < n_samples; i++) {
                        int32_t* old_ptr = vals + i * n_values_per_sample;
                        int32_t* new_ptr = new_vals + i * new_n_values;
                        int new_idx = 0;

                        // For Number=R, start from index 0 (ref); for Number=A, start from index 1 (first alt)
                        int start_idx = (id_len == BCF_VL_R) ? 0 : 1;

                        for (k = start_idx; k < n_allele_orig; k++) {
                            if (kbs_exists(keep_allele, k)) {
                                new_ptr[new_idx++] = old_ptr[k - start_idx];
                            }
                        }
                    }

                    bcf_update_format_int32(cfg.hdr_in, rec, key, new_vals, n_samples * new_n_values);
                    free(new_vals);
                } else if (!kbs_exists(cfg.warn_tag, hdr_id)) {
                    kbs_insert(cfg.warn_tag, hdr_id);
                    fprintf(stderr, "WARNING: got unexpected number of value from INFO/%s, so leaving it intact: (got %d, expected %d)\n", key, n_val, expected_values);
                }
            }
            if (vals)
                free(vals);
        } else if (field_type == BCF_HT_REAL) {
            int n_val = 0;
            float* vals = NULL;
            ret = bcf_get_format_float(cfg.hdr_in, rec, key, &vals, &n_val);
            if (ret > 0) {
                int n_values_per_sample = ret / n_samples;
                int expected_values = (id_len == BCF_VL_R) ? n_allele_orig : (n_allele_orig - 1);

                if (n_values_per_sample == expected_values) {
                    int new_n_values = (id_len == BCF_VL_R) ? n_keep : (n_keep - 1);
                    float* new_vals = (float*)malloc(n_samples * new_n_values * sizeof(float));

                    for (i = 0; i < n_samples; i++) {
                        float* old_ptr = vals + i * n_values_per_sample;
                        float* new_ptr = new_vals + i * new_n_values;
                        int new_idx = 0;

                        // For Number=R, start from index 0 (ref); for Number=A, start from index 1 (first alt)
                        int start_idx = (id_len == BCF_VL_R) ? 0 : 1;

                        for (k = start_idx; k < n_allele_orig; k++) {
                            if (kbs_exists(keep_allele, k)) {
                                new_ptr[new_idx++] = old_ptr[k - start_idx];
                            }
                        }
                    }

                    bcf_update_format_float(cfg.hdr_in, rec, key, new_vals, n_samples * new_n_values);
                    free(new_vals);
                } else if (!kbs_exists(cfg.warn_tag, hdr_id)) {
                    kbs_insert(cfg.warn_tag, hdr_id);
                    fprintf(stderr, "WARNING: got unexpected number of value from INFO/%s, so leaving it intact: (got %d, expected %d)\n", key, n_val, expected_values);
                }
            }
            if (vals)
                free(vals);
        } else if (field_type == BCF_HT_STR) {
            // String FORMAT fields with Number=R or Number=A are rare and complex to handle
            // For now, we'll just remove them
            if (!kbs_exists(cfg.warn_tag, hdr_id)) {
                kbs_insert(cfg.warn_tag, hdr_id);
                fprintf(stderr, "WARNING: Removing string multiallelic field INFO/%s as it will be inaccurate\n", key);
            }
            bcf_update_format_string(cfg.hdr_in, rec, key, NULL, 0);
        }
    }

    // Update genotypes
    int ngt = bcf_get_genotypes(cfg.hdr_in, rec, &ac, &n_ac);
    if (ngt > 0) {
        int max_ploidy = ngt / bcf_hdr_nsamples(cfg.hdr_in);
        int32_t* gt = ac; // Reuse the ac buffer for genotypes

        for (i = 0; i < bcf_hdr_nsamples(cfg.hdr_in); i++) {
            int32_t* ptr = gt + i * max_ploidy;
            for (j = 0; j < max_ploidy; j++) {
                // If the sample has fewer ploidy, the rest are bcf_int32_vector_end
                if (ptr[j] == bcf_int32_vector_end)
                    break;

                // Skip missing alleles
                if (bcf_gt_is_missing(ptr[j]))
                    continue;

                // Get the actual allele index (removing the bcf_gt_* flag)
                int allele = bcf_gt_allele(ptr[j]);

                // If this allele is being filtered out, set to missing
                if (allele >= n_allele_orig || !kbs_exists(keep_allele, allele)) {
#if 1
                    for (k = 0; k < max_ploidy; k++) {
                        ptr[k] = ptr[k] == bcf_int32_vector_end ? bcf_int32_vector_end : bcf_gt_missing;
                    }
                    break;
#else
                    ptr[j] = bcf_gt_missing;
                } else {
                    // Update the allele index to the new position
                    ptr[j] = bcf_gt_unphased(old2new[allele]);
#endif
                }
            }
        }

        bcf_update_genotypes(cfg.hdr_in, rec, gt, ngt);
    }
    // Update the record's alleles
    bcf_update_alleles(cfg.hdr_in, rec, (const char**)alleles, n_new);

    // Clean up
    free(alleles);
    free(old2new);
    kbs_destroy(keep_allele);
    if (ac)
        free(ac);
    if (af)
        free(af);
    if (ad)
        free(ad);

    cfg.n_sites_trimmed++;
    return rec;
}

/* Clean up. */
void destroy(void)
{
    kbs_destroy(cfg.warn_tag);
    fprintf(stderr, "Summary: Processed %" PRIu64 " sites, trimmed alleles at %" PRIu64 " sites\n",
        cfg.n_sites, cfg.n_sites_trimmed);
}
