/* @file main.c
 * @brief This file is a getoptslong example with Doxygen comments.
 *
 * @details
 * The details go here.
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

/* define the options */
const struct option long_options[] = {
    {"help",    no_argument,       0, 'h'},
    {"species", required_argument, 0, 's'},
    {"aligner", required_argument, 0, 'a'},
    {0, 0, 0, 0}
};

/* define the usage function */
const char *usage = "Usage: %s [options]\n"
		    "Options:\n"
		    "  -h, --help\t\tPrint this help message\n"
		    "  -s, --species\t\tSpecify the species\n"
		    "  -a, --aligner\t\tSpecify the aligner\n";
