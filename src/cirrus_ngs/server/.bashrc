# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

# User specific aliases and functions
alias getlogs="cd /shared/workspace/logs/ChipSeq/mus_test_proj/"
alias q="~/qstat"
alias c="clear"
alias pl="cd /shared/workspace/Pipelines/"
. /shared/workspace/software/software.conf
. ~/tempconfigure
