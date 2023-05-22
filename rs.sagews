from sage.misc.prandom import randint
from sage.coding.guruswami_sudan.interpolation import gs_interpolation_linalg
from sage.coding.guruswami_sudan.gs_decoder import roth_ruckenstein_root_finder

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

class GeneraliseReedSolomon:
    def __init__(self, F_q, k, x, y= None):
        if y == None:
            y= vector( [1]*len(x) )

        self.F_q= F_q
        self.k= k
        self.x= x
        self.y= y
        self.n= len(x)
        self.t= ((n - k)//2)

    def encode(self, m):
        C= codes.GeneralizedReedSolomonCode(self.x, self.k, self.y)

        # Generation de la matrice de Vandermonde du Code
        G = codes.encoders.GRSEvaluationVectorEncoder(C).generator_matrix()
#         G = codes.ReedSolomonCode(self.F_q, self.n, self.k).generator_matrix()

        return m * G

    def decode_bw(self, u):
        F_q= self.F_q
        t= self.t
        n= self.n
        M= []
        y= self.y

        for i in range(len(u)):
            u_i= u[i]
            x_i= self.x[i]

            v_exp1= [y[i] for j in range(n-t)]
            v_exp2= [(-1)* (x_i)^j for j in range(n-t)]
            v_i= [u_i * x_i^j for j in range(t+1)] + [v_i_exp1*v_i_exp2 for v_i_exp1,v_i_exp2 in zip(v_exp1, v_exp2)]

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

    def decode_sd(self, u, l, r):
        F_q= self.F_q
        t= self.t
        n= self.n
        _x= self.x

        R.<X,Y> = PolynomialRing(F_q, 'X,Y')
        pts= [ (F_q(_x[i]), F_q(u[i])) for i in range(n) ]

        Q= gs_interpolation_lee_osullivan(pts, l+1, (1, l), 1)
#         Q= gs_interpolation_linalg(pts, l+1, (1, 4), 1)
        Q= R(Q)

        roots= []
#         x = SR.var('x')
#         Q = Q.subs(X=x)
#         Q_y = Q.subs({'X': x})
#         print(Q)
#         Q_y= Q.change_ring(F_q['Y'])
#         print(type(Q_y))
#         print(SR(Q_y))

#         Q_y= Q_y.polynomial_expression()
#         print(Q_y)
#         print(L(Q_y))

#         roots = solve(Q_y == 0, Y)
#         roots = Q_y.roots()

#         Q_coefs = [Q.coefficient({y: i}) for i in range(Q.degree(y)+1)]
#         print(Q)
#         Q_y= sum(SR( Q_coefs[i] ) * y**i for i in range(len(Q_coefs)))
#         print(Q)
#         print(Q_coefs)
#         print(Q_y)
#         roots = Q_y.roots()
#         Q_coef = Q.terms('y')
#         Q_coef = Q_coef[::-1]

#         print(Q_y)
#         x = SR.var('x')
#         Q= Q.subs(x=x)
#         roots = solve(Q.subs(X=X), Y)
#         print(Q)
#         X = R.gens()[0]
#         Q_y = Q(X=X)
#         print(Q_y)
#         roots = Q_y.roots(multiplicities=False)

#         Q_y = Q.change_ring(R).subs('X', R.gen(0))
#         X, Y = SR.var('X Y')
#         Q_y = Q.change_ring(F_q)[0]
#         P = P.subs(x=2)

#         P.<X, Y> = PolynomialRing(F_q, 2)
#         Q= R(2*X + 3*X*Y + 3*Y^2)
#         print(Q)
#         roots = []

        def generate_polynomials(F_q, k):
            P = PolynomialRing(F_q, 'X')
            polynomials = []

            for deg in range(k):
                for coeffs in cartesian_product_iterator([F_q] * (deg + 1)):
                    poly = P(coeffs)
                    polynomials.append(poly)

            return polynomials

        polys= generate_polynomials(F_q, Q.degree(Y))
        for P in polys:
            if Q(X,P(X)) == 0:
                roots.append( P )

#         print(roots)
        return roots
#         Q_symbolic = Expression(SR, Q)
#         roots = solve(Q_symbolic, Y)
#         roots = solve(Q == 0, Y, algorithm='sympy')
#             print(Q_x)
#             Q_x_roots = Q_x.roots(multiplicities=False)
#         roots= Q.roots(Y)
#         roots = solve(Q, Y)

#             I = Ideal(Q_x)
#             variety = I.variety()
#             roots = [solution[y] for solution in variety]
#             print(variety)

#             roots.extend([(x, y) for y in Q_x_roots])

        print(roots)
#         roots = solve(P, P.variable(1))
#         x, y = P.parent().gens()
#         I = Ideal(P)
#         variety = I.variety()
#         roots = [solution[y] for solution in variety]
#         print(variety)
#         roots = P.roots(P.variable(1), field=False)
#         roots = solve(P, P.variable(1))
#         print(len(roots))
#         P= R(0)
#         for

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
t= ((n - k)//2)

F_q= GF(q)
canal= Canal(F_q, n, t)

x= list(F_q)[:n]
y= vector( [1]*len(x) )

# ReedSolomon : Berlekamp-Welch + Berlekamp Massey
# GeneraliseRS: Berlekamp-Welch

# Pour tester la version Generalize, décommenter la ligne suivante:
# y= list(F_q)[:n]
# y[0]= q-1

m= field_random_element_n(F_q, k)

print("ReedSolomon : Berlekamp-Welch + Berlekamp Massey")
print("ReedSolomon : Berlekamp-Welch", "\n")

print("Nombre de correction : ", t)
print("message clair : ", m, "\n")


RS= GeneraliseReedSolomon(F_q, k, x, y)
c_RS= RS.encode(m)
c_RS_sent= canal(c_RS)

x_RS_bw= RS.decode_bw(c_RS_sent)
x_RS_bm= RS.decode_bm(c_RS_sent)
x_RS_sd= RS.decode_sd(c_RS_sent, 5, t+2)

print("Message codé algo généré : ", c_RS)
print("Message codé algo envoyé : ", c_RS_sent)
print("Message algo Berl-Welch  : ", x_RS_bw)
print("Message algo Berl-Massey : ", x_RS_bm, "\n")
print("Message algo Sudan       : ", x_RS_sd, "\n")

# RS2 = codes.ReedSolomonCode(F_q, n, k)
RS2 = codes.GeneralizedReedSolomonCode(x, k, y)

c_RS2 = RS2.encode(m)
c_RS2_sent= canal(c_RS2)
x_RS2= RS2.decoder("BerlekampWelch").decode_to_message(c_RS2_sent)
# x_RS2= RS2.decoder("BerlekampWelch").decode_to_code(c_RS2_sent)

x_RS2= vector(x_RS2.list())
print("Message codé verif généré: ", c_RS2)
print("Message codé verif envoyé: ", c_RS2_sent)
print("Message codé verif reçu  : ", x_RS2, "\n")

print("--- Resultats ---")
print("Validité de Berlekamp-Welch  : ", x_RS_bw == m)
print("Validité de Berlekamp-Massey : ", x_RS_bm == m)
print("Validité de Sudan            : ", x_RS_sd == m)