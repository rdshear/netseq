version 1.0

import "netseq.wdl" as netseq_wf

workflow checker_netseq {
        meta {
        description: "Quick check for quay.io/rdshear/netseq workflow"
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    input {
        File inputFastQ 
        String refFasta
    }

    call netseq_wf.netseq {
        input:
            inputFastQ = inputFastQ,
            refFasta = refFasta
    }


}
