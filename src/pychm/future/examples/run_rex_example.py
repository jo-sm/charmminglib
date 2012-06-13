from pychm.future.scripts import rex_map

# the path to the log file
log = "/u/aokur/Ala10/fastrepd/repd.hel.100.rep10/replica.exch.gz"

# the directory where you want to write the unshuffled dcd files
out_dcd_dir = "~fpickard/hfrex_testing/asim_zip_test"

# a list, where each element is a path to an input dcd file, ORDERING IS IMPORTANT!
inp_dcd_fname = [
        "/u/aokur/Ala10/fastrepd/repd.hel.100.rep10/replica.dcd_0",
        "/u/aokur/Ala10/fastrepd/repd.hel.100.rep10/replica.dcd_1",
        "/u/aokur/Ala10/fastrepd/repd.hel.100.rep10/replica.dcd_2",
        "/u/aokur/Ala10/fastrepd/repd.hel.100.rep10/replica.dcd_3",
        "/u/aokur/Ala10/fastrepd/repd.hel.100.rep10/replica.dcd_4",
        "/u/aokur/Ala10/fastrepd/repd.hel.100.rep10/replica.dcd_5",
        "/u/aokur/Ala10/fastrepd/repd.hel.100.rep10/replica.dcd_6",
        "/u/aokur/Ala10/fastrepd/repd.hel.100.rep10/replica.dcd_7"
]

##########################################
# README #################################
##########################################
# the call signature of the rex_map function has changed!
# 1)
# some of the argument names have changed, the new arg names are below
# 2)
# a new optional argument `log_ftype` now exists.
# it defaults to 'auto', but also takes the following values:
# `None`: that is a literal `None` *not* the string "None"
#           `None means that this is an uncompressed logfile
# 'tar': archive with tar compression (untested)
# 'bz': archive with bz2 compression (untested)
# 'zip': archive with zip compression (untested)
# 'gz': archive with gzip compression (tested!)
# 'auto': guess a filetype based on "fname.endswith( ... )" not very robust...
# 3)
# the final argument is now a list instead of a `*list` as it was previously, just
# remove the extra asterisk from your call and it should work


rex_map(log_fname=log, out_dir=out_dcd_dir, log_ftype='auto', dcd_fnames=inp_dcd_fname)

