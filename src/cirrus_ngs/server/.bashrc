# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

# User specific aliases and functions
alias getlogs="cd /shared/workspace/logs/DNASeq/bwa_gatk/mus_test_proj/"
alias q="~/qstat"
alias c="clear"
alias pl="cd /shared/workspace/Pipelines/"
. /shared/workspace/Pipelines/config/software.conf

# added by Anaconda3 4.4.0 installer
export PATH="/shared/workspace/software/anaconda3/bin:$PATH"
