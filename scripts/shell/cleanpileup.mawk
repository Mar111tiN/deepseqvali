#!/bin/sh

### adjusted from mawktools to work on mpileup from one file and to ommit base qualities
### usage: samtools mpileup | cleanpileup
# removes from all reads the traces of base location and indel lengths (keeps IiDd as indel marker)


mawk '
{
    # print first 4 columns
    for (i=1;i++<4;) {
        printf("%s\t", $i);
    }
    read = $5;
    
    # remove position traces from all read fields
    gsub(/\^[^\t]|\$/,"",read);
    # 
    while (match(read,/[+-][0-9]+/)) {
        pre = substr(read,1,RSTART-2);
        indel = substr(read,RSTART,1); # + or -
        base = substr(read,RSTART-1,1); # base before the deletion
        l = substr(read,RSTART+1,RLENGTH-1);
        post = substr(read,RSTART+RLENGTH+l);
        if (indel == "-") {
            if (match(base,/[ACGT]/)) {
                base = "D";
            } else {
                base = "d";
            }
        } else {
            if (match(base,/[ACGT]/)) {
                base = "I";
            } else {
                base = "i";
            }            
        }

        read = pre base post;
    }     
# print all fields
    printf("%s\n",read);
}'
