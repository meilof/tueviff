# based on ...  Tomas Toft, CWI, TU/e SdH & BS
import math, operator, string, sys
from optparse import OptionParser
import viff.boostFP
viff.boostFP.install() 
from viff.runtimeFP import Runtime, create_runtime, gather_shares, Share
from viff.cdeferFP import setDebugging
from viff.comparisonFP import Toft07Runtime
from viff.config import load_config
from twisted.internet import reactor, defer
from twisted.python import log
import time
import viff.inlinecb



def load_tableau(filename):
    T = []
    comment_sign = "#"
    separator = "\t"
    for line in open("data/"+filename, "r"):
        # strip whitespace and comments
        tmp = string.split(line, comment_sign)[0]
        tmp = string.strip(tmp)
        if tmp and not tmp.startswith(comment_sign):
            # skip comments and empty lines.
            parsed_line = string.split(tmp, separator)
            T.append(map(int, parsed_line))
    T[-1].append(0)
    return T

    
@viff.inlinecb.inlinecb2
def solve(runtime, T):
    n = len(T[0]) - 1
    m = len(T) - 1
    for i in xrange(m+1):
        for j in xrange(n+1):
            T[i][j] = runtime.shareZp(T[i][j])
            
    yield [Share(runtime, runtime.Zp) for _ in xrange(n+m)]

    basis = [runtime.secret_index(i+n) for i in xrange(m)]
    cobasis = [runtime.secret_index(j) for j in xrange(n)]
    iter = 0
    prev_pivot = runtime.shareZp(1)
    
    while True:
        # test for termination
        print >>sys.stderr, "%d Termination?..." % iter
        min, p_col_index = runtime.int_minimal(T[-1][:-1])
        
        step = yield runtime.open(min < 0)
        if not step: break
        iter += 1
                
        # iterate
        print >>sys.stderr, "%d Determining pivot..." % iter
        p_col = [runtime.in_prod(T[i][:-1], p_col_index) for i in xrange(m+1)]
        constraints = [(T[i][-1] + (p_col[i] <= 0), p_col[i]) for i in xrange(m)]
        (_, pivot), p_row_index = runtime.frac_minimal(constraints)
        
        print >>sys.stderr, "%d Updating tableau..." % iter
        # swap basis entries
        minplus = map(operator.neg, p_row_index)+p_col_index
        row_col_ind = runtime.in_prod(basis+cobasis, minplus)
        h = runtime.scalar_mul(row_col_ind, minplus)
        for i in xrange(m):
            basis[i] -= h[i] 
        for j in xrange(n):
            cobasis[j] -= h[m+j]
            
        #  T[i,j] = T[i,j]*p/p' - (C[i]/p' - p_row_index[i])*(R[j] + p * p_col_index[j])
        prev_p_inv = ~prev_pivot
        p_col = runtime.scalar_mul(prev_p_inv, p_col)
        for i in xrange(m):
            p_col[i] -= p_row_index[i]
        p_row = [runtime.in_prod([T[i][j] for i in xrange(m)], p_row_index) for j in xrange(n+1)]
        delta_row = runtime.scalar_mul(prev_pivot, p_col_index)
        for j in xrange(n):
            p_row[j] += delta_row[j]
        T = runtime.gauss(T, pivot*prev_p_inv, p_col, p_row)
        prev_pivot = pivot
            
    print >>sys.stderr, "Termination..."
    max = yield gather_shares(map(runtime.open, [T[-1][-1],prev_pivot]))
    print >>sys.stderr, " max(f) = %d / %d" % (max[0].value, max[1].value)
    print >>sys.stderr, "Computing solution..."
    solution = runtime.secret_index_get(T[:-1][-1], basis, 0, n)
    print >>sys.stderr, "Computing dual solution..."
    dual_solution = runtime.secret_index_get(T[-1][:-1], cobasis, n, n+m)
    
    viff.inlinecb.returnValue(solution+dual_solution)

@viff.inlinecb.inlinecb2
def solve_lp(runtime, tableau_filename, certificate_filename):
    yield None
    
    print >>sys.stderr, "Initialization..."
    T = load_tableau(tableau_filename)
    n = len(T[0]) - 1
    m = len(T) - 1
    
    runtime.zp_secret_index(n+m)
    
    solop = yield gather_shares(map(runtime.open, solve(runtime, T)))
    
    print >>sys.stderr, "Writing solutions..."
    file = open("data/"+certificate_filename, "w")
    # #file.write('# tableau = \n' + options.tableau + '\n')
    # #file.write('# modulus = \n' + str(Zp.modulus) + '\n')
    # #file.write('# bit-length = \n' + str(runtime.options.bit_length) + '\n')
    # #file.write('# security param = \n' + str(runtime.options.security_parameter) + '\n')
    # #file.write('# threshold = \n' + str(runtime.threshold)+'\n')
    # #file.write('# common divisor = \n' + str(max[1].value) + '\n')
    file.write('# Solution = \n')
    for j in xrange(n):
        file.write(str(solop[j].value)+'\t')
    file.write('\n')
    file.write('# Dual Solution = \n')
    for i in xrange(m):
        file.write(str(solop[n+i].value)+'\t')
    file.write('\n')
    file.close()
    runtime.shutdown()


# Parse command line arguments.
parser = OptionParser()
parser.add_option("-f", "--tableau",
                  help="Filename for tableau.")
#parser.add_option("-t", "--threshold", dest="threshold", type="int",
#                  help="Threshold -- threshold should be below n/2")
parser.set_defaults(tableau="default.teau")
Runtime.add_options(parser)
options, args = parser.parse_args()

if len(args) == 0:
    parser.error("You must specify a config file.")
else:
    id, players = load_config(args[0])
if args[1]:
    certificate_filename = args[1]
else:
    print >>sys.stderr, "Please specify a filename to write certificate to:"
    certificate_filename = raw_input()

#log.startLogging(sys.stderr)

pre_runtime = create_runtime(id, players, options.threshold, options, Toft07Runtime)
pre_runtime.addCallback(solve_lp, options.tableau, certificate_filename)

reactor.run()
