# based on ...  Tomas Toft, CWI, TU/e SdH & BS
import operator, string, sys
from optparse import OptionParser
import viff.boost
viff.boost.install() 
from viff.field import GF
from viff.runtime import Runtime, create_runtime, gather_shares, Share
from viff.comparison import Toft07Runtime
from viff.config import load_config
from twisted.internet import reactor
from twisted.python import log
from gmpy import mpz

def load_tableau(filename):
    T = []
    comment_sign = "#"
    separator = "\t"
    for line in open(filename, "r"):
        # strip whitespace and comments
        tmp = string.split(line, comment_sign)[0]
        tmp = string.strip(tmp)
        if tmp and not tmp.startswith(comment_sign):
            # skip comments and empty lines.
            parsed_line = string.split(tmp, separator)
            T.append(map(int, parsed_line))
    T[-1].append(0)
    return T

def pow_list(a, x, n):
    if n==1:
        return [a]
    xs = []
    x2s = pow_list(a, x**2, n/2)
    for x2 in x2s:
        xs.append(x2)
        xs.append(x2*x)
    if n%2 == 1:
        xs.append(xs[-1]*x)
    return xs
    
def argmin(xs, gte_min):
    n = len(xs)
    if n == 1:
        return (xs[0], [1])
    if n == 2:
        c, min = gte_min(xs[0], xs[1])
        return (min, [1-c, c])
    x2s = []
    c2s = []
    while len(xs) > 1:       
        a1 = xs.pop()
        a0 = xs.pop()
        c, min = gte_min(a0, a1)
        x2s.append(min)
        c2s.append(c)
    if n%2 == 1:
        x2s.append(xs.pop())
    min, index2 = argmin(x2s, gte_min)
    index = []
    if n%2 == 1:
        index.append(index2[n/2])
    for i in xrange(n/2-1, -1, -1):
        b = index2[i]*c2s[i]
        index.append(index2[i] - b)
        index.append(b)
    return (min, index)

def int_minimal(xs):
    def int_gte_min(x0, x1):
        c = x0 >= x1
        return (c, c * (x1 - x0) + x0)
    return argmin(xs, int_gte_min)

class Protocol:

    def frac_minimal(self, xs):
        def frac_gte_min((n0, d0), (n1, d1)):
            c = self.runtime.in_prod([n0,d0],[d1,-n1]) >= 0
            h = self.runtime.scalar_mul(c, [n1-n0, d1-d0])
            return (c, (h[0]+n0, h[1]+d0))
        return argmin(xs, frac_gte_min)

    def __init__(self, runtime, tableau_filename, certificate_filename):
        print "Initialization..."
        self.runtime = runtime
        self.certificate_filename = certificate_filename
        self.T = load_tableau(tableau_filename)
        self.n = len(self.T[0]) - 1
        self.m = len(self.T) - 1
        k = runtime.options.security_parameter
        l = runtime.options.bit_length
        self.N = self.n + self.m
        if self.N%4 == 0: self.N += 1
        p = 1 + self.N*(1 + (self.N**2)%4 + 4*(2**(l+k+1)/self.N))
        while not mpz(p).is_prime():
            p += 4 * self.N
        self.Zp = GF(p)
        shareZp = lambda x:Share(self.runtime, self.Zp, self.Zp(x))
        for i in xrange(self.m+1):
            for j in xrange(self.n+1):
                self.T[i][j] = shareZp(self.T[i][j])
        w = self.Zp(2)**((p-1)/self.N)  # lucky choice 2
        self.w_powers = [self.Zp(1)]
        for _ in xrange(self.N-1):
            self.w_powers.append(self.w_powers[-1] * w)
        self.basis = [Share(self.runtime, self.Zp, self.w_powers[-(i+self.n)]) for i in xrange(self.m)]
        self.cobasis = [Share(self.runtime, self.Zp, self.w_powers[-j]) for j in xrange(self.n)]
        self.iter = 0
        self.prev_pivot = shareZp(1)
        self.iterate()
        
    def iterate(self):
        print "%d Termination?..." % self.iter
    	min, index = int_minimal(self.T[-1][:-1])
        step = self.runtime.open(min < 0)
        step.addCallback(self.test_termination, index)
        
    def test_termination(self, step, index):
        if step:
            self.iter += 1
            self.prev_pivot = self.iteration(index)
            self.iterate()
        else:
            self.terminate()   

    def iteration(self, p_col_index):
        print "%d Determining pivot..." % self.iter
        p_col = [self.runtime.in_prod(self.T[i][:-1], p_col_index) for i in xrange(self.m+1)]
        constraints = [(self.T[i][-1] + (p_col[i] <= 0), p_col[i]) for i in xrange(self.m)]
        (_, pivot), p_row_index = self.frac_minimal(constraints)
			
        print "%d Updating tableau..." % self.iter
        # swap basis entries
        minplus = map(operator.neg, p_row_index)+p_col_index
        row_col_ind = self.runtime.in_prod(self.basis+self.cobasis, minplus)
        h = self.runtime.scalar_mul(row_col_ind, minplus)
        for i in xrange(self.m):
            self.basis[i] -= h[i] 
        for j in xrange(self.n):
            self.cobasis[j] -= h[self.m+j]
            
        #  T[i,j] = T[i,j]*p/p' - (C[i]/p' - p_row_index[i])*(R[j] + p * p_col_index[j])
        prev_p_inv = ~self.prev_pivot
        p_col = self.runtime.scalar_mul(prev_p_inv, p_col)
        for i in xrange(self.m):
            p_col[i] -= p_row_index[i]
        p_row = [self.runtime.in_prod([self.T[i][j] for i in xrange(self.m)], p_row_index) for j in xrange(self.n+1)]
        delta_row = self.runtime.scalar_mul(self.prev_pivot, p_col_index)
        for j in xrange(self.n):
            p_row[j] += delta_row[j]
        self.T = self.runtime.gauss(self.T, pivot*prev_p_inv, p_col, p_row)
        return pivot
            
    def terminate(self):
        print "Termination..."
        f = gather_shares(map(self.runtime.open, [self.T[-1][-1],self.prev_pivot]))
        f.addCallback(self.print_max)
        
    def print_max(self, max):
        print " max(f) = %d / %d" % (max[0].value, max[1].value)
        print "Computing solution..."
        solution = [0]*self.n
        for i in xrange(self.m):
            x_powers = pow_list(self.T[i][-1]*~self.Zp(self.N), self.basis[i], self.N)
            for j in xrange(self.n):
                coefs = [self.w_powers[(j*k)%self.N] for k in xrange(self.N)]
                solution[j] += self.runtime.lin_comb(coefs,x_powers)
        print "Computing dual solution..."
        dual_solution = [0]*self.m
        for j in xrange(self.n):
            x_powers = pow_list(self.T[-1][j]*~self.Zp(self.N), self.cobasis[j], self.N)
            for i in xrange(self.m):
                coefs = [self.w_powers[((self.n+i)*k)%self.N] for k in xrange(self.N)]
                dual_solution[i] += self.runtime.lin_comb(coefs,x_powers)
        solutions = gather_shares(map(self.runtime.open, solution+dual_solution))
        solutions.addCallback(self.write_solutions, max[1].value)

    def write_solutions(self, solutions, cd):
        file = open("data\\"+self.certificate_filename, "w")
        file.write('# tableau = \n' + options.tableau + '\n')
        file.write('# modulus = \n' + str(self.Zp.modulus) + '\n')
        file.write('# bit-length = \n' + str(self.runtime.options.bit_length) + '\n')
        file.write('# security param = \n' + str(self.runtime.options.security_parameter) + '\n')
        file.write('# threshold = \n' + str(self.runtime.threshold)+'\n')
        file.write('# common divisor = \n' + str(cd) + '\n')
        file.write('# Solution = \n')
        for j in xrange(self.n):
            file.write(str(solutions[j].value)+'\t')
        file.write('\n')
        file.write('# Dual Solution = \n')
        for i in xrange(self.m):
            file.write(str(solutions[self.n+i].value)+'\t')
        file.write('\n')
        file.close()
        self.runtime.shutdown()

# Parse command line arguments.
parser = OptionParser()
parser.add_option("-f", "--tableau",
                  help="Filename for tableau.")
parser.add_option("-t", "--threshold", dest="threshold", type="int",
                  help="Threshold -- threshold should be below n/2")
parser.set_defaults(tableau="default.teau", threshold = 1)
Runtime.add_options(parser)
options, args = parser.parse_args()

if len(args) == 0:
    parser.error("You must specify a config file.")
else:
    id, players = load_config(args[0])
if args[1]:
    certificate_filename = args[1]
else:
    print "Please specify a filename to write certificate to:"
    certificate_filename = raw_input()

log.startLogging(sys.stdout)

pre_runtime = create_runtime(id, players, options.threshold, options, Toft07Runtime)
pre_runtime.addCallback(Protocol, options.tableau, certificate_filename)

reactor.run()
