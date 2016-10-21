#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>

#include "tree.h"
#include "treestats.h"

struct treekernel_options {
    double decay_factor;
    double gauss_factor;
    double sst_control;
    int nLTT;
    int normalize;
    int ladderize;
    scaling scale_branches;
    FILE *tree1_file;
    FILE *tree2_file;
};

struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"decay-factor", required_argument, 0, 'l'},
    {"gauss-factor", required_argument, 0, 'g'},
    {"sst-control", required_argument, 0, 's'},
    {"nLTT", no_argument, 0, 'c'},
    {"normalize", no_argument, 0, 'n'},
    {"ladderize", no_argument, 0, 'd'},
    {"scale-branches", required_argument, 0, 'b'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: treekernel [options] [tree1] [tree2]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
    fprintf(stderr, "  -l, --decay-factor        penalty for large matches (default 0.2)\n");
    fprintf(stderr, "  -g, --gauss-factor        variance for Gaussian RBF (default 2)\n");
    fprintf(stderr, "  -c, --nLTT                multiply by nLTT (default no)\n");
    fprintf(stderr, "  -s, --sst-control         0 for subtree kernel, 1 for subset-tree kernel (default 1)\n");
    fprintf(stderr, "  -n, --normalize           divide output by square-root of product of\n");
    fprintf(stderr, "                            kernels of trees with themselves\n");
    fprintf(stderr, "  -d, --ladderize           ladderize trees before computing kernel\n");
    fprintf(stderr, "  -b, --scale-branches      type of branch scaling to apply (mean/median/max/none, default none)\n");
}

struct treekernel_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct treekernel_options opts = {
        .decay_factor = 0.2,
        .gauss_factor = 2,
        .sst_control = 1.0,
        .nLTT = 0,
        .normalize = 0,
        .ladderize = 0,
        .scale_branches = NONE,
        .tree1_file = stdin,
        .tree2_file = stdin
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "hl:g:s:cndb:", long_options, &i);
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 'l':
                opts.decay_factor = atof(optarg);
                break;
            case 'g':
                opts.gauss_factor = atof(optarg);
                break;
            case 's':
                opts.sst_control = atof(optarg);
                break;
            case 'c':
                opts.nLTT = 1;
                break;
            case 'n':
                opts.normalize = 1;
                break;
            case 'd':
                opts.ladderize = 1;
                break;
            case 'b':
                if (strcmp(optarg, "mean") == 0) {
                    opts.scale_branches = MEAN;
                }
                else if (strcmp(optarg, "median") == 0) {
                    opts.scale_branches = MEDIAN;
                }
                else if (strcmp(optarg, "max") == 0) {
                    opts.scale_branches = MAX;
                }
                break;
            case '?':
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
        }
    }

    // TODO: safety
    if (optind < argc) {
        opts.tree1_file = fopen(argv[optind++], "r");
    }
    if (optind < argc) {
        opts.tree2_file = fopen(argv[optind++], "r");
    }

    return opts;
}

int main (int argc, char **argv)
{
    struct treekernel_options opts = get_options(argc, argv);
    double knum, kdenom = 1;

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_t *t1 = parse_newick(opts.tree1_file);
    igraph_t *t2 = parse_newick(opts.tree2_file);


    if (opts.ladderize) {
        ladderize(t1);
        ladderize(t2);
    }

    scale_branches(t1, opts.scale_branches);
    scale_branches(t2, opts.scale_branches);

    if (opts.normalize) {
        kdenom = sqrt(kernel(t1, t1, opts.decay_factor, opts.gauss_factor, opts.sst_control)) *
                 sqrt(kernel(t2, t2, opts.decay_factor, opts.gauss_factor, opts.sst_control));
    }

    knum = kernel(t1, t2, opts.decay_factor, opts.gauss_factor, opts.sst_control);
    if (opts.nLTT) {
        knum *= 1.0 - nLTT(t1, t2);
    }

    printf("%f\n", knum / kdenom);

    igraph_destroy(t1);
    igraph_destroy(t2);

    if (opts.tree1_file != stdin) {
        fclose(opts.tree1_file);
    }
    if (opts.tree2_file != stdin) {
        fclose(opts.tree2_file);
    }
    return EXIT_SUCCESS;
}
