#!/bin/sh

# 1 CORE

### cleans samtools mpileup output
# for deletions, start is raised by one while assuming that coords had been lowered in bed file to capture deletion

mawk ' # cat - <(echo "STOP\n") $SNP |   ###### if the SNP markers are necessary, they have to be read in here 
## $5 ~ /[.,]*[ACTG]+[.,]*/ | mawk    # if you only want to output mutant positions (complete the quotes!!)


NR == 1 { 
    minVAF="'${1-0}'";
    print("minVAF:", minVAF) >> "/dev/stderr";
    ###### QUERY ##############
    # get the letters to look for
    len = split("Aa,Gg,Cc,Tt,Dd,Ii",Letters,",")
    ###### HEADER ################
    for (l=0;l++<len;) {
        LETPAT[l] = "[" Letters[l] "]"
    }
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Chr","Start","Ref", "Alt", "Depth", "TR2", "VAF");
    ####### HAVE TO ADJUST ###########
    # loop through the letters
}
# special case: no coverage
$4 == 0 {
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $5, $4,$4, -1);
    next;
}

{   ######### LINES #############
    # loop through the letters
    for (l = 0;l++< len;) {
        COUNT[l] = gsub(LETPAT[l], "", $5);
    }
  ######### OUTPUT #############
  # loop through the letters
    maxcount = 0;
    base = "";
    vaf = 0;
    ref = $3;
    start = $2;
    # first line extra for pretty
    for (l=0;l++<len;) {
        if (COUNT[l] > maxcount) {
            maxcount = COUNT[l];
            base = substr(Letters[l],1,1);
            vaf =maxcount/$4;
        } 
    }
    if (match(base, "D")) {
        # for deletions readjust the base to "-"
        # readjust the ref to "?" because here base is the one preceeding the deletion
        # reup the start coords that had been decremented for bed coords
        base = "-";
        ref = "?";
        start++;
    } 

    if (vaf== 0) {
        base = ref;
    }
    if (vaf >= minVAF) {
        # I want chr in chromosomes
        if ($1 !~ "chr") {
            $1 = "chr" $1;
        }
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, start, ref, base,$4, maxcount, vaf);
    }
}'