
from collections import namedtuple
import gmpy2 as gy 
from random import randint
import matplotlib.image as image
import matplotlib.pyplot as plt
import hashlib
import math


Point = namedtuple("Point", "x y")
O = 'Origin'

p = 4001 #9 bit
a = 2
b = 3

def modp(a,p):
	return a%p

def eq_modp(a,b,p):
	return modp(a-b,p)==0

def inverse_modp(x,p):
	for y in range(p):
		if eq_modp(x * y, 1, p):
			return y
	return None

def valid(P,p,a,b):
    if P == O:
        return True
    else:
        return (eq_modp((P.y**2),(P.x**3 + a*P.x + b),p) and (0 <= P.x < p) and (0 <= P.y < p))

def ec_inv(P,p):
    if P == O:
        return P
    return Point(P.x, modp(-P.y,p))


def ec_add(P,Q,a,b,p):
    # if not (valid(P,p,a,b) and valid(Q,p,a,b)):
    #     raise ValueError("Invalid inputs")
    if P == O:
        result = Q
    elif Q == O:
        result = P
    elif Q == ec_inv(P,p):
        result = O
    else:
        if P == Q:
            dydx = (3 * P.x**2 + a) * inverse_modp(2 * P.y,p)
        else:
            dydx = modp((Q.y - P.y),p) * inverse_modp((Q.x - P.x),p)
        v=modp((P.y - dydx*P.x),p)
        x = modp((dydx**2 - P.x - Q.x),p)
        y = modp((-dydx * x - v),p)
        result = Point(x, y)
    # assert valid(result,p,a,b)
    return result

def ec_multiply(P,k,a,b,p):
	res=O
	for i in range(k):
		if res==O:
			res=P
		else:
			res=ec_add(res,P,a,b,p)
	return res

def order_E(a,b,p):
	points=[]
	points.append(O)
	for i in range(p):
		for j in range(p):
			if eq_modp(j**2,(i**3)+a*i+b,p):
				points.append(Point(i,j))
	return (len(points),points)


def choose_n(N):
	k=1
	while(k<=N):
		k=gy.next_prime(k)
		if N%k==0:
			n=k
	return n
def cofactor(N,n):
	return N//n

def choose_base_point(h,points,a,b,p):
	for i in points:
		G=ec_multiply(i,h,a,b,p)
		if(G!=O):
			break
	return G 

		 
(N,points)=order_E(a,b,p)
#print(points)
n=choose_n(N)
h=cofactor(N,n)
P=choose_base_point(h,points,a,b,p) #base_point
print("Order of Elliptic Curve = ",N)
print("Order of Subgroup = ",n)
print("Cofactor of Subgroup = ",h)
print("Base Point of Subgroup = ",P)


#user S=Alice
#user R=Bob
def dec2bin(a,size):
	b=bin(a).replace("0b", "")
	return '0'*(size-len(b))+b 

def gen_keys(p,a,b,G,n,h):
	private_key=randint(1,n-1)
	public_key=ec_multiply(G,private_key,a,b,p)
	return (private_key,public_key)


def sha_256_en(P,Zsr,p,a,b):
	input_image=image.imread('input.jpg')
	h=int(hashlib.sha256(input_image.tobytes()).hexdigest(),16)
	#print('31')
	e=ec_multiply(P,h%p,a,b,p)
	#print('32')
	he=int(hashlib.sha256(str(e.x^e.y).encode()).hexdigest(),16)
	#print('33')
	Ke=ec_add(e,Zsr,a,b,p)
	#print('34')
	H=int(hashlib.sha256(str( Ke.x ^ Ke.y ).encode()).hexdigest(),16)
	#print('35')
	return(Ke,H,he)


def chaotic_map(H,l):
	s0=2
	s1=4
	Hb=bin(H).replace("0b","")
	Hd=[]
	for i in range(0,len(Hb),8):
		Hd.append(int(Hb[i:i+7],2))
	sth0=s0+(sum(Hd[0:16])/256)
	sth1=s1+(sum(Hd[16:32])/256)
	ch_map=[]
	ch_map.append(int(sth0))
	ch_map.append(int(sth1))
	for i in range(2,l):
		ch_map.append(math.ceil(abs(0.8-(10*(math.sin(3.14*ch_map[i-1])**2))+(0.01*1*abs(1-(2*ch_map[i-2]))))))
	return ch_map


def image_processing(size,a,b,p):
	input_image=image.imread('input.jpg')
	plt.figure()
	plt.imshow(input_image,cmap="gray",vmin=0,vmax=255)
	plt.title("Input Image")
	l=[]
	pr_image=[]
	for i in range(len(input_image)):
		for j in range(len(input_image[0])):
			l.append(dec2bin(input_image[i][j],size))
			pr_image.append(input_image[i][j])
	#print(len(pr_image),"pr")
	im_points=[]
	for i in range(len(pr_image)):
		x=pr_image[i]
		x=x*14
		while(int(math.sqrt(((x**3)+a*x+b)%p)+0.5)**2 != ((x**3)+a*x+b)%p):
			x=x+1
		y=((x**3+a*x+b)%p)**0.5
		im_points.append(Point(x,y))
	#print(len(im_points),"im_points")

	return im_points

def encry(P,a,b,p):
	Vs,Zs=gen_keys(p,a,b,P,n,h)
	Vr,Zr=gen_keys(p,a,b,P,n,h)
	Zsr=ec_multiply(P,Vs*Vr,a,b,p)
	#print("1")
	points=image_processing(8,a,b,p)
	#print("2")
	Ke,H,he=sha_256_en(P,Zsr,p,a,b)
	#print("3")

	L=[]
	for i in range(len(points)):
		L.append(ec_add(Ke,points[i],a,b,p))
	S1=[]
	for i in L:
		S1.append(int(i.x))
		S1.append(int(i.y))
	#print(len(S1),"s1")
	PR=chaotic_map(H,len(S1))
	C=[]
	for i in range(len(S1)):
		C.append(S1[i]^PR[i])
	return C,Ke,H


Ct,Ket,Ht=encry(P,a,b,p)
print("Cphiper_Image:")
print(Ct)
#print(len(Ct),"en")


def decry(P,Ke_r,H_r,C,a,b,p):
	PR=chaotic_map(H_r,len(C))
	R=[]
	for i in range(len(C)):
		R.append(C[i]^PR[i])
	#print(len(R),"dec_R")
	DR=[]
	for j in range(0,len(R)-1,2):
		DR.append(Point(R[j],R[j+1]))
	#print(len(DR),"dec_DR")
	PT=[]
	for k in range(len(DR)):
		PT.append(ec_add(DR[k],Point(Ke_r.x,-Ke_r.y),a,b,p))
	#print(len(PT),"dec_PT")
	MS=[]
	for i in PT:
		MS.append(int(i.x/14))
	#print(len(MS),"dec_MS")
	Output_Image=[]
	for i in range(0,len(MS),256):
		Output_Image.append(MS[i:i+255])
	plt.figure() 
	plt.imshow(Output_Image,cmap="gray",vmin=0,vmax=255)
	plt.title("decypted Image")
	plt.show()
decry(P,Ket,Ht,Ct,a,b,p)





