# Fibanocci numbers module 

class Fib(object):

    def __init__(self, n):
        self.n = n

    def fib(self): # write Fibanocci series up to n 
        a, b = 0,1 
        while b < self.n:
            print b
            a, b = b, a+b 
      
    def fib2(self): # return Fibanocci series up to n
        result = [] 
        a, b = 0,1 
        while b < self.n:
            result.append(b)
            a,b = b, a+b 
        return result 

#if __name__ == "__main__":
#    import sys
#    fib(int(sys.argv[1]))
