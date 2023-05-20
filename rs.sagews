from sage.misc.prandom import randint

import random

def field_random_element_distinct(K, s):
    x= s
    while x == s:
        x= K.random_element()
    return x

def field_random_element_n(K, n):
    x= []
    while n:
        x.append( K.random_element() )
        n-= 1
    return vector(x)

class Canal:
    def __new__(cls, K, n, t):
        def q_canal(x):
            y= []
            if len(x) != n:
                print("Erreur de taille de n")
                return

            can= random.sample([i for i in range(n)], t)
            for i in range(0, n):
                if i in can:
                    y.append( field_random_element_distinct(K, x[i]) )
                else:
                    y.append( x[i] )
            return vector(y)

        return q_canal

#     def send_canal(x):
#         return self(x)

class ReedSolomon:
    def __init__(self, F_q, k, x):
        self.F_q= F_q
        self.k= k
        self.x= x
        self.n= len(x)
        self.t= ((n - k)//2)

    def encode(self, m):
        C= codes.GeneralizedReedSolomonCode(self.x, self.k)

        # Generation de la matrice de Vandermonde du Code
        G = codes.encoders.GRSEvaluationVectorEncoder(C).generator_matrix()
#         G = codes.ReedSolomonCode(self.F_q, self.n, self.k).generator_matrix()
        return m * G

    def decode_bw(self, u):
        F_q= self.F_q
        t= self.t
        n= self.n
        M= []

        for i in range(len(u)):
            u_i= u[i]
            x_i= self.x[i]
            v_i= [u_i * x_i^j for j in range(t+1)] + [(-1)* (x_i)^j for j in range(n-t)]
            M.append(v_i)

        R= PolynomialRing(F_q, 'x')

        v_i= [0]*t+[1]+[0]*(n-t)
        M.append(v_i)

        M= matrix(F_q, M)
        b= column_matrix(vector([0]*(n)+[1])).change_ring(F_q)
        V= M.augment(b)
        V= V.rref()

        A= R( V[:t+1, -1].list() )
        B= R( V[t+1:, -1].list() )

        F, rem = B.quo_rem(A)

        if rem != 0:
            print("Erreur de décodage avec Berlekamp-Welch")
            return -1

        return vector( F.list() )

    def decode_bm(self, u):
        F_q= self.F_q
        t= self.t
        n= self.n
        _x= self.x

        R= PolynomialRing(F_q, names='x')
        x = R.gen()

        pts = [(_x[i], u[i]) for i in range(n)]
        U = R.lagrange_polynomial(pts)

        P= R(1)

        for x_i in _x:
            P*= (x - x_i)

#         error_poly = -U
#         print(error_poly)

        def ext_euclide(U):
            C= [ R(1), R(0), R(0) ]
            A= [ R(0), 1   , R(0) ]
            B= [ P   , U   , R(0) ]

            i= 0
            while B[1] != 0:
                Q, B[2]= B[0].quo_rem(B[1])
                A[2] = A[0] - Q * A[1]
                C[2] = C[0] - Q * C[1]


                i+= 1
                if A[2].degree() <= t and B[2].degree() < n - t:
                    return A[2], B[2], C[2]


                # Shifting our vars

                A[0]= A[1]
                A[1]= A[2]

                B[0]= B[1]
                B[1]= B[2]

                C[0]= C[1]
                C[1]= C[2]


            return 0

        A, B, C= ext_euclide(U)
        F, rem = B.quo_rem(A)

        if rem != 0:
            print("Erreur de décodage avec Berlekamp-Massey")
            return -1

        return vector( F.list() )

#         y= []
#         for i in range( len(u) ):
#             x_i= self.x[i]
#             if A(x_i) == 0:
#                 y.append(F( x_i ))
#             else:
#                 y.append(u[i])

#         return vector(y)

q= 11
n= 7
k= 3
# k= n - t*2
t= ((n - k)//2)

F_q= GF(q)
canal= Canal(F_q, n, t)

# x= [F_q(i) for i in range(n)]
x= list(F_q)[:n]
m= field_random_element_n(F_q, k)

print("Nombre de correction : ", t)
print("message clair : ", m, "\n")


RS= ReedSolomon(F_q, k, x)
c_RS= RS.encode(m)
c_RS_sent= canal(c_RS)

x_RS_bw= RS.decode_bw(c_RS_sent)
x_RS_bm= RS.decode_bm(c_RS_sent)

print("Message codé algo généré : ", c_RS)
print("Message codé algo envoyé : ", c_RS_sent)
print("Message algo Berl-Welch  : ", x_RS_bw)
print("Message algo Berl-Massey : ", x_RS_bm, "\n")

# RS2 = codes.ReedSolomonCode(F_q, n, k)
RS2 = codes.GeneralizedReedSolomonCode(x, k)

c_RS2 = RS2.encode(m)
c_RS2_sent= canal(c_RS2)
x_RS2= RS2.decoder("BerlekampWelch").decode_to_message(c_RS2_sent)
# x_RS2= RS2.decoder("BerlekampWelch").decode_to_code(c_RS2_sent)

x_RS2= vector(x_RS2.list())
print("Message codé verif généré: ", c_RS2)
print("Message codé verif envoyé: ", c_RS2_sent)
print("Message codé verif reçu  : ", x_RS2, "\n")

print("--- Resultats ---")
print("Validité de Berlekamp-Massey : ", x_RS_bm == m)









