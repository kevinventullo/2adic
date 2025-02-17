from sys import argv
from functools import cmp_to_key

class TwoAdic:
    def __init__(self):
        pass
    
    
    def toInt(self, x):
        rv = 0
        for i, c in enumerate(reversed(x)):
            if (c != '0' and c != '1'):
                raise Exception('Improperly formatted argument: ', x, '. Detected character ', c)
            if c == '1':
                rv += 2**i
        return int(rv)
    
    # Returns a string of the first n digits of x in base 2
    # If n is not provided, returns a string of all digits, stopping at leading 0's or 1's
    def toBitString(self,x, n = None):
        rv = ''
        if n is None:
            if x == 0:
                return '0'
            if x == -1:
                return '1'
            while(x!=0 and x != -1):
                rv = str(x%2) + rv
                x = x // 2
        else:
            for _ in range(n):
                rv = str(x%2) + rv
                x = x // 2
        return rv
    
    # Compute x! modulo n
    def modFac(self, x, n):
        rv = 1
        for i in range(2,x+1):
            rv  = (rv*i)%n
        return int(rv)
    
    # # Compute x^k modulo n; repeated squaring could speed this up, or even other methods in the blogpost.
    # def modExp(self, x, k, n):
    #     rv = 1
    #     for _ in range(k):
    #         rv = (rv * x)%n
    #     return int(rv)
    
    def TwoAdicInv(self, x, n = None):
        #print('Computing inv of ', x, ' with ', n, ' digits specified')
        if(x%2) == 0:
            raise Exception('Only handles two-adic integers for now')
        if n is None:
            n = len(self.toBitString(x))
        # Simple linear. Can replace w/ Newton-Raphsom to get lg(n) time. 
        rv = x
        for i in range(1,n+1):
            if ((rv*x)%(2**i) != 1):
                rv = rv + (2**(i-1))
        #print('Returning ', int(rv))
        return int(rv)
    
    def TwoAdicDiv(self, x, y, n = None):
        #print('Computing div of ', x, ' and ', y, ' with ', n, ' digits specified')
        
        if (y == 0):
            raise Exception('Divided by zero!')
        if n is None:
            n = len(self.toBitString(x))
        while(x % 2 == 0 and y % 2 == 0):
            x = x // 2
            y = y // 2
        if (y % 2 == 0):
            raise Exception('Only handles two-adic integers for now.')
        #print('No really computing div of ', x, ' and ', y, ' with ', n, ' digits specified')
        rv = int((x*self.TwoAdicInv(y, n))%(2**n))
        #print('Returning ', rv)
        return rv
    
    def TwoAdicLog(self,x, n = None):
        if(x%4 != 1):
            raise Exception('Arguments to log must be 1 mod 4')
        if n is None:
            # Determine how many digits to compute based on size of x
            n = len(self.toBitString(x))
        x = x - 1
        rv = 0
        # Since x is 1 mod 4, the nth term of the log Taylor series is always divisible by at least 2**n
        # Thus, only need to compute up to the n-1th term. 
        for i in range(1,n):
            rv += (-1)**(i+1) * self.TwoAdicDiv(x**i, i, n)
            rv = rv % (2**n)
            #print('Intermediate value: ', rv)
        if rv < 0:
            rv = rv + 2**n
        
        return int(rv)
    
    def TwoAdicExp(self,x, n = None):
        if(x%4 != 0):
            raise Exception('Arguments to exp must be 0 mod 4')
        if n is None:
            # Determine how many digits to compute based on size of x
            n = len(self.toBitString(x))
        rv = 1
        # Since x is 0 mod 4, the nth term of the exp Taylor series is always divisible by at least 2**n
        # Thus, only need to compute up to the n-1th term. 
        for i in range(1,n):
            rv += self.TwoAdicDiv(x**i, self.modFac(i, 2**n), n)
            rv = rv % (2**n)
        if rv < 0:
            rv = rv + 2**n
        return int(rv)
    

    def TwoAdicStringLog(self, x, n = None, abbrev = False):
        if n is None:
            n = len(x)
        if abbrev:
            x = x + '01'
            n = n + 2
        x_int = self.toInt(x)
        rv_int = self.TwoAdicLog(x_int, n)
        rv = self.toBitString(rv_int, n)
        if abbrev:
            rv = rv[:-2]
        
        return rv
    
    def TwoAdicStringExp(self, x, n = None, abbrev = False):
        if n is None:
            n = len(x)
        if abbrev:
            x = x + '00'
            n = n + 2
        x_int = self.toInt(x)
        rv_int = self.TwoAdicExp(x_int, n)
        rv = self.toBitString(rv_int, n)
        if abbrev:
            rv = rv[:-2]
        return rv


TA = TwoAdic()
#print(TA.toBitString(5))
m = 6
n = 7
if len(argv) > 1:
    m = int(argv[1])
if len(argv) > 2:
    n = int(argv[2])
else:
    n = m+1
for i in range(m,n):
    print('Log table for the kth digit of log; k = ', i)
    tmplist = []
    for j in range(2**i):
        jstr = TA.toBitString(j, i)
        jlog = TA.TwoAdicStringLog(jstr, abbrev=True)
        if jstr[0] == '0' and jlog[0] == '1':
            x = ''
            if (jlog[-1] == '0' or jlog[-3] == '1'):
                x = '0'
            else:
                x = '1'
            y = ''
            if (jlog[-4] == '1' and jlog[-2] == '1'):
                y = '1'
            else:
                y = '0'
            tmplist.append((jstr[1:], jlog[0:], x, y))
    def compare(item1, item2):
        if item1[1] < item2[1]:
            return -1
        if item1[1] > item2[1]:
            return 1
        return 0
    
    #tmplist = sorted(tmplist, key=cmp_to_key(compare))
    for elt in tmplist:
        #print(jstr[1:], jlog[0:], x, y)
        print(elt[0], elt[1], elt[2], elt[3])
    
        

#print(TA.toBitString(TA.TwoAdicLog(5, 5), 5))