

import os.path
from collections import deque

from pychm.future.tools import paragraphs, _myexpandpath, _myopenzip_context
from pychm.future.io import open_dcd


class LogError(Exception):
    """Exception raised when exchange log parsing goes wrong."""
    pass


class LogEntry(object):
    """
    """
    __slots__ = ["exchange", "step", "repeat", "body", "index_array", "temp_array"]

    def __init__(self):
        self.exchange = None
        self.step = None
        self.repeat = None
        self.body = None
        self.index_array = None
        self.temp_array = None

    def _parse_body(self):
        self.index_array = tuple(int(line.split()[0]) for line in self.body )
        self.temp_array = tuple(float(line.split()[1]) for line in self.body )

    def __repr__(self):
        return "%s(%r, %r, %r)" % (self.__class__.__name__, self.exchange, self.step, self.repeat)

    def __len__(self):
        return len(self.body)


class ExchangeLog(object):
    """
    """
    def __init__(self, fname, validate=True, ftype=None):
        self.text = self._init_text(fname, ftype)
        self._validate_text()
        self.header = self.text.popleft()
        self.log = self.parse_log()
        if validate:
            self._validate_log()

    @property
    def temp_array(self):
        try:
            return self.log[0].temp_array
        except IndexError:
            return None

    def _init_text(self, fname, ftype):
        """open input file, normalize the text, dump it to a deque"""
        with _myopenzip_context(fname, ftype) as inp_fp:
            iterator = ( line.lower() for line in inp_fp )
            iterator = ( line.strip() for line in iterator )
            iterator = ( line for line in iterator if line )
            return deque(iterator)

    def _validate_text(self):
        """check if the header exists, and check if exchange delimeters exist"""
        if not self.text[0].startswith('#'):
            raise LogError('Bad formatting, missing log file header, should be one line and start with "#".')
        if not self.text[1].startswith('# exchange'):
            raise LogError('Bad formatting, missing exchange delimeter, should be one line and start with "# Exchange".')

    def _validate_log(self):
        """checks if all log entries are the same length, and have the same
        temp_array the way the temp_array is checked is pretty lazy. this might
        have to change to account for floating point rounding. running this
        method in the constructor is the default, but can be shut off by
        setting the keyword `validate=False`"""
        try:
            length = len(self.log[0])
        except IndexError:
            raise LogError("Can't determine length of first log element.")
        try:
            temp = self.log[0].temp_array
        except AttributeError:
            raise LogError("First log element has no temp_array attribute")
        for i, entry in enumerate(self.log):
            if length != len(entry):
                raise LogError("Element %d in log does not have the correct length of %d" % (i, length))
            if temp != entry.temp_array:
                raise LogError("Element %d in log does not have the correct temperature array of %r" % (i, temp))

    def parse_log(self):
        """read through self.text, instantiate LogEntry objects and return a
        list of them"""
        def _parse_entry(iterable):
            tmp = LogEntry()
            # parse values from exchange header
            header = iterable[0].split()
            tmp.exchange = int(header[2][:-1])
            tmp.step = int(header[4][:-1])
            tmp.repeat = int(header[6])
            tmp.body = iterable[1:]
            tmp._parse_body()
            return tmp
        tmp = []
        for entry in paragraphs(self.text, "# exchange"):
            tmp.append(_parse_entry(entry))
        return tmp

    def __iter__(self):
        return iter(self.log)

    def __len__(self):
        return len(self.log)

def get_nsavc(*fnames):
    if not fnames:
        raise ValueError("must specify input dcd files")
    nsavc_array = []
    for fname in fnames:
        with open_dcd(fname, mode='r') as fp:
            nsavc_array.append(fp.nsavc)
    nsavc0 = nsavc_array[0]
    for nsavc in nsavc_array:
        if nsavc != nsavc0:
            raise ValueError('bad value for nsavc, first dcd has nsavc = %d but dcd file "%s" has value of %d' % (nsavc, dcd_fnames[i], nsav))
    return nsavc0

###############################################################################
# work time ###################################################################
###############################################################################
def rex_map(log_fname, out_dir, *dcd_fnames):
    # validate inputs
    if not isinstance(log_fname, basestring):
        raise TypeError("invalid log_fname")
    if not isinstance(out_dir, basestring):
        raise TypeError("invalid out_dir")
    if not dcd_fnames:
        raise ValueError("must specify input dcd files")
    # parse rexlog
    rexlog = ExchangeLog(log_fname)
    # parse nsavc value
    nsavc = get_nsavc(*dcd_fnames)
    # make nsavc is compatible with step values from log
    for entry in rexlog:
        if entry.step % nsavc != 0:
            raise ValueError("bad value for nsavc, is not multiple of step %d in exchangelog entry %d" % (entry.step, entry.exchange))
    # create output fnames
    out_dcd_fnames = []
    for i, temp in enumerate(rexlog.temp_array):
        fname = _myexpandpath(out_dir) + os.path.sep + 'rex_map_%02d_%6.2f.dcd' % (i, temp)
        out_dcd_fnames.append(fname)
    # write headers
    with open_dcd(dcd_fnames[0], mode='r') as fp:
        base_header = fp.export_header()
    for fname in out_dcd_fnames:
        with open_dcd(fname, mode='wb') as fp:
            tmp_header = {}
            tmp_header.update(base_header)
            tmp = ["\x02\x00\x00\x00"]
            tmp.append("%-80s" % "* This DCD was generated by rex_map")
            tmp.append("%-80s" % "* This DCD index is %d, and temp is %d" % (i,rexlog.temp_array[i]))
            tmp_header['title'] = ''.join(tmp)
            fp.import_header(tmp_header)
            fp.write_header()
    # write frames
    inp_dcd = [ open_dcd(fname, mode='rb') for fname in dcd_fnames ]
    out_dcd = [ open_dcd(fname, mode='ab') for fname in out_dcd_fnames ]
    ## ignore entries where repeat != 1
    iterator = ( entry for entry in rexlog if entry.repeat == 1 )
    current_step = 0 # how many steps have we written?
    for entry in iterator:
#        print "exchange number: %d" % entry.exchange
        while current_step < entry.step:
#            print "  step number: %d" % (current_step + nsavc)
            ###################
            # do work #########
            ###################
            for i, out in enumerate(out_dcd):
#                print "    reading from: %s" % dcd_fnames[entry.index_array[i]-1]
#                print "    writing to: %s" % out.name
                out.write(inp_dcd[entry.index_array[i]-1].read_frame())
            ##################
            # done with work #
            ##################
            current_step += nsavc

###############################################################################
# inputs ######################################################################
###############################################################################

if __name__ == '__main__':
    inp_fname = "/u/tim/projects/hfreq_rex/test1/rex_standard_repeat.exch"
    inp_dcd_fname = ["/u/tim/projects/hfreq_rex/test1/rex_standard_repeat.dcd_0",
            "/u/tim/projects/hfreq_rex/test1/rex_standard_repeat.dcd_1",
            "/u/tim/projects/hfreq_rex/test1/rex_standard_repeat.dcd_2",
            "/u/tim/projects/hfreq_rex/test1/rex_standard_repeat.dcd_3"]
    out_dcd_dir = "/u/fpickard/hrex_testing/tim_test/"
    rex_map(inp_fname, out_dcd_dir, *inp_dcd_fname)
