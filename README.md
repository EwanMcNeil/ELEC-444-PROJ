

CONTENTS

[5](#br11)

[6](#br11)

[7](#br11)

[8](#br11)

[Original-M3](#br11)[ ](#br11). . . . . . . . . . . . . . 11

[Voxel](#br11)[ ](#br11)[Selection-M1](#br11)[ ](#br11). . . . . . . . . . 11

[Voxel](#br11)[ ](#br11)[Selection-M2](#br11)[ ](#br11). . . . . . . . . . 11

[Voxel](#br11)[ ](#br11)[Selection-M3](#br11)[ ](#br11). . . . . . . . . . 11

[BlockWise-M1](#br12)[ ](#br12). . . . . . . . . . . . . 12

[BlockWise-M2](#br12)[ ](#br12). . . . . . . . . . . . . 12

[BlockWise-M3](#br12)[ ](#br12). . . . . . . . . . . . . 12

[Classic](#br13)[ ](#br13)[Computation](#br13)[ ](#br13)[Time-M1](#br13)[ ](#br13). . . . 13

[Voxel](#br13)[ ](#br13)[Selection](#br13)[ ](#br13)[Computation](#br13)[ ](#br13)[Time](#br13)[ ](#br13)[-M1](#br13)[ ](#br13)13

[I](#br2)

[Abstract](#br2)

2

2

[II](#br2)

[III](#br3)

[Introduction](#br2)

[Methods](#br3)

[9](#br12)

[10](#br12)

[11](#br12)

[12](#br13)

[13](#br13)

3

3

[III-A](#br3)

[III-B](#br3)

[Classical](#br3)[ ](#br3)[Implementation](#br3)[ ](#br3).

[Automatic](#br3)[ ](#br3)[Tuning](#br3)[ ](#br3)[of](#br3)[ ](#br3)[the](#br3)

[Smoothing](#br3)[ ](#br3)[parameter](#br3)[ ](#br3). . .

[Voxel](#br3)[ ](#br3)[Selection](#br3)[ ](#br3). . . . . . .

[Block-Wise](#br4)[ ](#br4). . . . . . . . .

[threading](#br5)[ ](#br5). . . . . . . . . .

3

4

4

5

[III-C](#br4)

[III-D](#br4)

[III-E](#br5)

[IV](#br5)

[Implementations](#br5)

5

5

5

5

6

[IV-A](#br5)

[Classical](#br5)[ ](#br5). . . . . . . . . . .

[IV-A.1](#br5)

[IV-A.2](#br5)

[IV-A.3](#br6)

[getNewValue](#br5)[ ](#br5). .

[weight](#br5)[ ](#br5). . . . . .

[Intensities](#br5)[ ](#br5). . . .

[IV-B](#br6)

[Classical](#br6)[ ](#br6)[with](#br6)[ ](#br6)[Voxel](#br6)[ ](#br6)[Selec-](#br6)

[tion](#br6)[ ](#br6). . . . . . . . . . . . . .

6

[IV-B.1](#br7)

[IV-B.2](#br7)

[IV-B.3](#br8)

[GetNewValueVoxel](#br7)[ ](#br7)7

[ZVoxelSelection](#br7)[ ](#br7).

[WVoxelSelection](#br8)

7

8

8

9

[IV-C](#br8)

[IV-D](#br9)

[Blockwise](#br8)[ ](#br8). . . . . . . . . .

[Multi-Processor](#br9)[ ](#br9). . . . . . .

[V](#br10)

[Results](#br10)

[V-A](#br10)

[V-B](#br11)

10

[Classical](#br10)[ ](#br10). . . . . . . . . . . 10

[Classical](#br11)[ ](#br11)[with](#br11)[ ](#br11)[Voxel](#br11)[ ](#br11)[Selec-](#br11)

[tion](#br11)[ ](#br11). . . . . . . . . . . . . . 11

[Blockwise](#br12)[ ](#br12). . . . . . . . . . 12

[Multiprocessing](#br12)[ ](#br12). . . . . . . 12

[V-C](#br12)

[V-D](#br12)

[VI](#br13)

[Problems](#br13)[ ](#br13)[faced](#br13)

[VI-A](#br13)

[VI-B](#br13)

13

[Size](#br13)[ ](#br13)[of](#br13)[ ](#br13)[input](#br13)[ ](#br13)[image](#br13)[ ](#br13). . . . 13

[Selection](#br13)[ ](#br13)[of](#br13)[ ](#br13)[correct](#br13)[ ](#br13)[nor-](#br13)

[malizing](#br13)[ ](#br13)[constant](#br13)[ ](#br13)[Z](#br13)[ ](#br13). . . . 13

[Python](#br13)[ ](#br13)[negative](#br13)[ ](#br13)[indexing](#br13)[ ](#br13). 14

[Blockwise](#br14)[ ](#br14). . . . . . . . . . 14

[VI-C](#br14)

[VI-D](#br14)

[VII](#br14)[ ](#br14)[Conclusion](#br14)

[VIII](#br14)[ ](#br14)[References](#br14)

14

14

LIST OF FIGURES

[1](#br10)

[2](#br10)

[3](#br10)

[4](#br11)

[Ground](#br10)[ ](#br10). . . . . . . . . . . . . . . . . 10

[Noisy](#br10)[ ](#br10). . . . . . . . . . . . . . . . . . 10

[Original-M1](#br10)[ ](#br10). . . . . . . . . . . . . . 10

[Original-M2](#br11)[ ](#br11). . . . . . . . . . . . . . 11





NL-Mean Denoising Implementation

Lorenzo Pes Ewan McNeil

Gina Cody School of Electrical and Computer Engineering - 2020

Course: ELEC 444 - Medical Image Processing

I. ABSTRACT

possible; Non Local means simulates this concept

by averaging parts of the image with other parts

that are similar to it. This differs from traditional

”local” denoising algorithms by searching through

the whole volume of the image for similar patches

to average with in comparison to averaging with

neighboring pixels.

A common problem in medical image process-

ing of Magnetic Resonance Images (MRI) is that

of removing noise components present at the out-

put image of the machine. This process is named

”denoising” and one of its crucial aspect is that

of being able to remove noise components from

the image without degrading relevant information

such as edges. One of the most common and

effective methods for denoising MRI images is

called Non-Local Means. The objective of this

algorithm is that of exploiting self similarities

within an image to remove noise components.

Given that MRI images are commonly in a 3D

format the computational complexity of the al-

gorithm becomes quite extensive. For that reason

four incremental improvements were proposed to

reduce the complexity of the algorithm while keep-

ing relevant information intact. Therefore, the main

purpose of this report is to implement the solutions

proposed by Coupe (2008) to cope with the high

computational burden of the algorithm. These im-

provements are ﬁrst the Automatic Tuning of the

Smoothing parameter, Voxel Selection,Blockwise,

and Multi threading. Of the four improvements

it was found that multi-threading and Block wise

allowed for large impacts on the time complexity,

while improvements like Voxel Selection increased

the edge clarity.

By using this method, the end result compared

to local denoising methods is an increase in ﬁl-

tering of noise and an increase in edge clarity.

However, with these ﬁltering improvements comes

a huge computational complexity. With ”true” Non

Local means, every pixel is compared to every

other pixel in the image and with a large enough

image, this will take a long period of time. This

paper will go into the theory of four improvements

upon the Non Local means algorithm proposed by

Coupe[1], and then the Python implementation of

these improvements created for this project.

The ﬁrst improvement is the automatic tuning

of the smoothing parameter, which increases the

clarity of the image and ensures that the neigh-

borhood size doesn’t affect the outcome. The sec-

ond improvement proposed by Coupe[1] is Voxel

selection wherein only voxels that will have an

impact on the averaging are calculated as opposed

to comparing all the voxels that are within the

current voxels search volume thus reducing the

complexity of the algorithm.

The Blockwise improvement is an abstraction

on-top of Non Local means where voxels are

II. INTRODUCTION

The Non Local means algorithm is an important grouped into overlapping blocks and the algorithm

algorithm that is used in various applications from is run on the blocks, thus allowing the Non Local

denoising in medical imaging modalities to in- means algorithm to run in fewer cycles. The voxel

painting (removing objects) from photos. Concep- values are then calculated based on all the out-

tually, a very good way of denoising an image is to comes of blocks that overlap it.

take multiple shots of the same image and average

Lastly, the fourth improvement upon the Non

them, thus achieving closer to the ground truth. Local means ﬁlter is recognizing that one iteration

In many applications, however, this simply isn’t of the algorithm doesn’t depend on the outcome





of a previous iteration. Therefore, the algorithm

can be run in parallel on subsections of the whole

input image, allowing for a drastic reduction in

time complexity. These improvements demonstrate

ways in which the complexity of the algorithm can

be reduced, thereby making the algorithm more

efﬁcient and thus more effective.

After a detailed discussion of the theory behind

the implementation of Non Local means and its

improvements, this paper will illustrate and explain

the implementations of each improvement iteration

in Python and then compare and contrast the

results obtained with the results from Coupe[1].

Lastly, the paper will discuss the problems faced

and the resolutions contained within each iteration

of the improvements.

X

w(x , x ) = 1

i

(3)

j

x∈Ω3

From equation (1) it can be deduced that each

each voxel at xi is compared with all the re-

maining voxels xj in the entire volume of the

image. Clearly, this is impractical because the time

complexity would be huge. For that reasons the

number of voxels to which xi is compared are

only limited to a speciﬁc volume called ”search

volume” (Vi) of size (2M+1)3. Therefore, for each

voxel x in V the weighted Euclidean distance

j

j

between the neighbor pixels N of voxel x and the

i

i

neighbor pixels N of voxel x of size (2d + 1)3

j

j

is calculated. Consequently, the weight function is

calculated as follow:

III. METHODS

1

||

−

||2

u(N

) u(N )

The methods section will be divided up ac-

cording to the improvements that the paper pro-

poses. Most of these improvements build upon

one another and thus there are overlaps in method

implementations or slight changes to speciﬁc pre-

vious methods to full-ﬁll these improvements. The

paper will start by explaining the methods used

in the Classical implementation of the Non-Local

Means De-noising algorithm and then examine the

changes to the methods when implementing the

Voxel Selection, Block wise and Multi threading

improvements to the time complexity of the algo-

rithm proposed by the paper.

i

j

w(x , x ) =

i

e

2

(4)

h

j

Zi

where Zi is a normalization constant that ensure

that (2) remains valid. On the other hand, the

parameter h2 is the smoothing parameters that

controls the decay of the exponential function.

The smoothing parameter ensures that the restored

value of voxel xi is the weighted average of

the voxels xi with most similar neighbor. It is

important to note that the smoothing parameter

is dependent from the neighbor size. In the next

section a method for making h independent of the

neighbor size.

A. Classical Implementation

B. Automatic Tuning of the Smoothing parame-

In the classical implementation of the NL-

mean algorithm the restored intensity NL(u)(xi)

of voxel xi is the weighted average of all voxels

intensities in the image u. This weighted average

is deﬁne as follow:

ter

As shown in equation (4) the smoothing param-

eter h is dependent from the neighbor size Ni. The

main objective of obtaining a smoothing parameter

that is independent of the nighbor size is that

of establishing the relationship between the noise

X

w(x , x )u(x )

j

(1)

NL(u)(xi) =

i

j

standard deviation σ2, the neighbor size

Ni

the constant β. Given that white Gaussian noise

and

x∈Ω3

where u(x ) is the intensity of voxel x and Ω3 is a good approximation for the real noise, the

j

j

is the entire volume of the 3D image. Speciﬁcally, standard deviation can be estimated by pseudo-

residuals using the following formula:

w(x , x ) is the weight function used to quantify

i

j

the similarity between voxels x and x under the

constraints that:

i

j

r

X

6

7

1

6

i =

(u(xi)) −

u(xj)

(5)

w(x , x ) ∈ [0, 1]

(2)

i

j

x∈P

i





Where P are the 6-neighbor at voxel x . The D. Block-Wise

i

i

standard deviation can be found as a least square

estimator as follows:

The Block-Wise improvement upon the Classi-

cal Implementation allows the algorithm to yield

a very similar output while drastically reducing

the spacial complexity. It does this by abstracting

the input image voxels into overlapping blocks,

doing non-local means on these blocks and then

extracting the values of the voxels by averaging

the all the block values that overlap it. In order

for this functionality to work the size of the block

is declared as ”a” and the stepping distance or

the distance between two centers of blocks is

N. In order to ensure that the blocks do overlap

then the equality 2a > N must hold true. After

setting these two parameters the Non-Local means

weighting equation is applied to the blocks in a

similar fashion as in the Classic implementation

the equation for this is deﬁned as follows:

X

1

σ2 =

2

i

(6)

|Ω3|

i∈Ω3

Furthermore, in order to make the ﬁlter indepen-

dent of Ni and to reduce complexity, instead of

using a weighted Euclidean distance, a classical

Euclidean Distance will be used. In conclusion,

equation (4) will take the following form:

2

1

||u(N )−u(N )||

i

j

w(x , x ) =

i

e

2βσ2|N

|

(7)

i

j

Zi

Where only parameter β becomes an hyper-

parameter to be tuned. It is important to note that

in our implementation this weight function will be

used for the classical implementation rather than

(4).

X

NL(u)(Bik) =

w(B , B )u(B ) (9)

ik j j

x∈Ω3

Where the only key difference between Blockwise

and the classic is that in the Blockwise weighting

summation we are dealing with blocks that contain

many voxels instead of comparing a single voxel

to another single voxel.The equation for the weight

follows again a very similar format to the classical

implementation but again deals with block com-

parisons and not voxel ones. The slightly differing

weight equation is show below.

C. Voxel Selection

Voxel Selection is the second method proposed

to reduce computational complexity. The main pur-

pose of this improvement is to select only the most

relevant voxels in the search volume beforehand. In

other words, only the voxels that have the highest

weight and therefore are the most similar will be

used in the computation of the Euclidean Distance.

In order to achieve the previous objective, the a

priori selection of voxels xj will be based one the

mean variance of u(N ) and u(N ). Finally, the

2

1

||u(B )−u(B )||

ik

j

w(B , B ) =

j

e

2βσ |N

2

|

(10)

i

ik

Zi

Which follows the same logic as the previous

weighting equation in classic and so will not be

explained again, however with Blockwise there is

a third step that needs to occur which is unpacking

the voxel values from the blocks that overlap them,

mathematically this is show as follows.

i

j

weight function will be calculated as follows:





||u(Ni)−u(Nj)||2

2|

w(x , x ) =

i

2βσ Ni

j

1

j

(8)

1 e

|

if µi

<

u(Ni) < 1 andσ

<

i

V ar(u(Ni)) <

1

Z

u(N

)

µ

V ar(u(N ))

σ2

i

j

1



0 otherwise

X

1

NL(u)(xi) =

Ai(P)

(11)

Ai

In the equation above, u(N ) and V ar(u(N ))

P∈A

i

i

i

are the mean and variance of the neighbor of voxel

Which deﬁnes the outputs voxel value as the

xi respectively. In the code implementation of this average of all the block value vectors that overlap

part, it will be seen that there were some issues in that voxel. The number of overlaps can vary due

the actual implementation in particular to ﬁnd the to the voxels position relative to the center of the

normalizing constant Zi.

blocks and can also vary due to the initial Stepping





value of N and block size of A that where deﬁned

\1) getNewValue: The main purpose of this

beforehand in the Blockwise implementation. With function is that of generating the search volume for

an increased block size there will be more blocks a given M. The function takes 3 input parameters:

that overlap on a single voxel and thus more

vectors for that voxel to average and get its new

ﬁltered value.

intuple: voxel

x

.

i

•

• indata: the input image u.

E. threading

• Z: normalizing constant.

Non local means is a non iterative algorithm

meaning that future calculations don’t depend on

previous calculations. For example the calculation important to note that Z is a constant for classical

for a voxel in the middle of the image could implementation since we can ﬁnd its value a priori:

The code for this function is shown below. It is

theoretically done before the calculation for the

voxel at the top left of the image where the

algorithm normally begins. This property of the

algorithms allows sections of the image to be run

in parallel due to the fact that the outputs of the al-

gorithm don’t impact other output. The paper ”An

optimized Blockwise Non-local means denoising

ﬁlter for 3-D magnetic resonance images” suggests

to ”divide the volume of the whole image into sub

volumes” and then run the algorithm on these sub

volumes. In the implementation proposed here in

this report the volumes that are deﬁned are slices of

the image that the algorithm then runs on the multi-

threading and pipe lining capabilities that Python

libraries provide.

Z = ((2M) + 1)3

def getNewValue(intuple, indata ,Z):

total = 0;

suma = 0;

global M

\# this deﬁnes the search volume

for x in range( intuple [0]−M,intuple[0]+M+1):

for y in range( intuple [1]−M,intuple[1]+M+1):

for z in range( intuple [2]−M,intuple[2]+M+1):

w = weight(intuple ,( x,y,z) , indata ,Z)

total = total + w indata[x,y,z]

\*

suma = suma + w

return total

In the code above, the for loop is generating the

IV. IMPLEMENTATIONS

set of tuple x to which the intuple (x ) generated

All the implementations of the algorithm were

made using Python language, Numpy library to

treat the image as a matrix and MatplotLib library

to plot obtained results.

j

i

in the main function is being compared to. The

variable w is holding the returned value of weight

between x and x . The variable total is fundamen-

i

j

tally implementing equation (1). It is important to

note that the for loop is taking into consideration

that the indata was padded by size M before

being processed. The working principle of function

weight will be explained in the next subsection.

Moreover, the padding function is shown below.

A. Classical

In the classical implementation equation (7)

is implemented for the weight of equation (1).

Multiple function were deﬁned to implement the

previously mentioned functions. The list of func-

tions is as follow:

def padding(image,m):

padded = np.pad(image, (( m, m), (m, m), (m,

m)),mode=’constant’, constant values =0)

return padded

• getNewValue (intuple,indata,Z)

• weight (tupleI,tupleJ,data,Z)

• Intensities (a,d, tuplein)

\2) weight: The main purpose of this function is

that of generating the weight expressed in equation

An explanations for each of the functions is (7). The function takes 4 input parameters:

given in the incoming subsections.

• tupleI: voxel xi.





• tupleJ: voxel x .

\+ data small [index7,index2,index4]

\+ data small [index7,index3,index4]

\+ data small [index7,index1,index5]

\+ data small [index7,index1,index6]

j

• data: image u.

• Z: normalizing constant.

The implementation of the function is as fol-

lows:

except IndexError:

total = total +0

def weight( tupleI , tupleJ , data ,Z):

epsilon [x,row,column] = math.pow(

math.sqrt (6/7)

(

global s

b = 1

N = 27

\*

data small [x,row,column] − (

(1/6) total ) ),2 )

\*

output = 0

return ( 1/( data small . size ) ) np.sum(epsilon)

\*

uNJ = np.asarray (

On the other hand, b represents β explained in

section 3B. In order to implement the Euclidean

distance linalg.norm function belonging to Numpy

library was used. Furthermore, the variable expo-

nential holds the exponential part of equation (4).

It is important to note that uNj and uNi holds the

neighbor intensities of voxel x and x respectively,

neighboorhoodIntensities ( data ,1, tupleJ ) )

uNI = np.asarray (

neighboorhoodIntensities ( data ,1, tupleI ) )

dist = np. linalg .norm(uNJ−uNI)

div = dist /(2 b s N)

\* \* \*

exponential = math.exp((−1) div)

\*

output = (1/Z) exponential

\*

i

j

return output

returned by the function Intensities(data,1, tupleI).

In the code above, s represents equation (6)

found using the following implementation:

\3) Intensities: The function is implemented as

follow.

def sigma(data ) :

data small = data

epsilon = np. zeros ( data small .shape)

def Intensities (a,d, tuplein ):

return [[[ a[x][y][z] if not(x==tuplein[0] and

y==tuplein[1] and z==tuplein [2]) and x >= 0

for x in range( data small .shape[0]) :

for row in range( data small .shape[1]) :

for column in range( data small .shape[2]) :

index1=row

< > <

and x a.shape[0] and y = 0 and y

a.shape[1] and z = 0 and z a.shape[2] else 0

\>

<

for x in range( tuplein [0]−d, tuplein [0]+1+d)]

for y in range( tuplein [1]−d,

tuplein [1]+1+d)]

index2=row+1

index3=row−1

index4=column

index5=column+1

for z in range( tuplein [2]−d,

tuplein [2]+1+d)]

index6=column−1

index7=x

index8=x+1

index9=x−1

total =0

try :

if index6<0 or index3<0 or

index9<0 :

total =

data small [index8,index1,index4]

\+ data small [index7,index2,index4]

\+ data small [index7,index1,index5]

else :

This function return matrix a[x,y,z] containing all

the intensities of the neighbor pixels of tuplein.

The variable d represents the neighbor size. The

comparisons statements ensures that if the voxel

under consideration is close to the edge, then

a value of 0 is selected for that neighbor. This

statements become obsolete if the image is padded,

but the padding was implemented afterwards. They

could be removed in later versions of the imple-

mentation.

B. Classical with Voxel Selection

In this section, equation (8) will be imple-

mented. It is important to note that for this section

total =

data small [index8,index1,index4]

\+ data small [index9,index1,index4]





the previously deﬁned function Intensities will

remain the same. The only function that will be

changed is that used to calculate the weight and to

get the new value. Furthermore, another function

to compute the normalization constant Zi will be

implemented since it can not be know a priori

\# this deﬁnes the search volume

for x in range( intuple [0]−M,intuple[0]+M+1):

for y in range( intuple [1]−M,intuple[1]+M+1):

for z in range( intuple [2]−M,intuple[2]+M+1):

w = WVoxelSelection(intuple ,( x,y,z) , indata

,Z)

what is the number of voxels that will be selected.

The implementation of the previous created a huge

overhead that will results in computational time

that are actually longer that the classic implemen-

tation. In further version of the implementation this

function and the location where it is called must

be rethought to avoid the increase in computational

time. Therefore the following new function intro-

duced are the following:

total = total + w indata[x,y,z]

\*

suma = suma + w

return total

\2) ZVoxelSelection: The purpose of this func-

tion is that of calculating the number of selected

voxel in order to compute the value of the nor-

malizing constant Z. The function is very similar

to WVoxelSelection with the only difference that

in the ﬁrst only a count of the number of selected

voxel is returned. The implementation is as fol-

lows:

• GetNewValueVoxel(intuple, indata)

• ZVoxelSelection(intuple,(x,y,z),indata)

• WVoxelSelection(intuple,(x,y,z),indata,Z)

def ZVoxelSelection( tupleI , tupleJ , data ) :

An explanations for each of the functions is

given in the incoming subsections.

u1 = 0.95

s1 = 0.5

\1) GetNewValueVoxel: The main purpose of

this function is the same as that of section A.1.

That is, to create the search volume for each voxel

xi generated in main. The main difference from

the previous is that we need to calculated the

normalization constant beforehand. As it can be

seen in the implementation in the next page. A for

loop that runs across the search volume is used to

ﬁnd the total number of Z based on the number

of selected voxels. Given that this for loop is run

z = 0

uNJ = np.asarray ( Intensities ( data ,1, tupleJ ) )

uNI = np.asarray ( Intensities ( data ,1, tupleI ) )

uNI mean var=uNI

uNJ mean var=uNJ

uNI mean var[1,1,1] = data [ tupleI ]

uNJ mean var[1,1,1] = data [ tupleJ ]

every time a new voxel xi is fed to the function

the computational time increase. This is deﬁnitely

mean i, var i = N var mean(uNI mean var)

mean j, var j = N var mean(uNJ mean var)

a wrong implementation and it must be revisited

in further sections.

if ( (mean i==0 and var i==0) or (mean j==0 and

var j==0) ):

div mean = math.inf

div var = math.inf

def GetNewValueVoxel(intuple, indata ):

else :

total = 0;

suma = 0;

Z = 0

div mean = mean i/mean j

div var = var i / var j

global M

if ( ( (div mean > u1)

and (div mean < 1/u1) )

and ( ( div var > s1)

and ( div var < 1/s1) ) ):

z = z +1

for x in range( intuple [0]−M,intuple[0]+M+1):

for y in range( intuple [1]−M,intuple[1]+M+1):

for z in range( intuple [2]−M,intuple[2]+M+1):

Z = Z

+ZVoxelSelection( intuple ,( x,y,z) , indata )

else :





z = z + 0

return z

if ( ( (div mean > u1) and (div mean < 1/u1) ) and

( ( div var > s1) and ( div var < 1/s1) ) ):

dist = np. linalg .norm(uNJ−uNI)

The variables s1 and u1 are corresponding to the

mean and variance limits deﬁned in equation 8.

The mean and variance for a given neighborhood

is found using the Numpy library in the following

function.

div = dist /(2 b s N)

\* \* \*

exponential = math.exp((−1) div)

\*

output = (1/Zv) exponential

\*

return output

else :

return 0

def N var mean(N):

return np.mean(N),np.var(N)

C. Blockwise

Where N is the neighborhood of either voxel xi or

xj previously found using Intensities function.

The code implementation of Blockwise is split

into two functions the ﬁrst function getNewVal-

ueBlock() implements the Blockwise version of

getNewValue() seen before in the previous im-

provements of the code. It implements equation

Eight seen in the Blockwise explanation where all

the blocks in the current blocks search volume

are compared against and the values multiplied

by their weight and summed with a normalization

factor. The Python implementation of this function

is seen as follows:

\3) WVoxelSelection: This function is used to

determine the weight of the voxels that are selected

by the bound of mean and variance speciﬁed

in variable s and u . The implementation is as

1

1

follows:

def WeightVoxelSelection( tupleI , tupleJ , data ,Zv):

global s

global M

b = 1

N = 27

def getNewValueBlock(intuple, indata ,M,N):

u1 = 0.95

s1 = 0.5

total = 0;

suma = 0;

count = 0

#blockZ = (N 2)

blockZ = 8

\# this deﬁnes the search and only compares on

blocks with the step of N

for x in range( intuple [0]−M,intuple[0]+M+1,N):

for y in

output = 0

3

\* \*\*

uNJ = np.asarray (

neighboorhoodIntensities ( data ,1, tupleJ ) )

uNI = np.asarray (

neighboorhoodIntensities ( data ,1, tupleI ) )

range( intuple [1]−M,intuple[1]+M+1,N):

for z in

uNI mean var=uNI

uNJ mean var=uNJ

range( intuple [2]−M,intuple[2]+M+1,N):

count = count + 1

uNI mean var[1,1,1] = data [ tupleI ]

uNJ mean var[1,1,1] = data [ tupleJ ]

w = weight(intuple

,( x,y,z)

, indata ,blockZ)

total = total + w indata[x,y,z]

suma = suma + w

mean i, var i = N var mean(uNI mean var)

mean j, var j = N var mean(uNJ mean var)

\*

return total

if ( (mean i==0 and var i==0) or (mean j==0 and

var j==0) ):

The method that is calling this function is re-

ferred to as blockwise() wherein the blocks are

inferred to by their center voxel and the stepping

range, therefore there is no need for an overhead

of a data structure to hold their values. The ﬁrst

div mean = math.inf

div var = math.inf

else :

div mean = mean i/mean j

div var = var i / var j





part of the method passes every block is passed speedup in the time it took for the algorithms to

through the get new value block where their values complete. In order to implement this the Python

are saved in an array that is the size of the output library ”multiprocessing” this library was chosen

array. The values of each block are mapped to instead of the ”Threading” library due to the

the location of the center voxel deﬁning them. fact that multiprocessing creates sub processes to

Therefor this matrix allows for an easy way to run instead of threads this, allows for individual

get the voxel values in the second part of the processes on separate processors. The threading

blockwise() function by simply inferring all the library on the other hand may create threads that

blocks that over lap the voxel given the A value. are still running on the same processor but that

The implementation of this code is seen as follows: is multi threading with context switching for the

implementation the goal was to have dedicated

def blockwise(M,padded,data):

N = 2

processors running these processes rather than

threads which may be sharing a processor. With

the Multiprocessor goal in mind the classical Non-

local means functions have to be altered ever so

slightly in order for them to be run properly in

a multiprocessing pool. Shown below are the two

implementations of the parallel functions for the

classic and the voxel selection.

A = 1

output = np. zeros ( data .shape)

centers = []

\## this gets all the center values

for x in range(M,padded.shape[0]−M,N):

for y in range(M,padded.shape[1]−M,N):

for z in range(M,padded.shape[2]−M,N):

output[x−M,y−M,z−M] =

getNewValueBlock((x,y,z)

,padded,M,N)

def parallelize classic (x):

centers .append((x−M,y−M,z−M))

output thread = np. zeros ( data .shape)

for y in range(M,padded.shape[1]−M):

for z in range(M,padded.shape[2]−M):

output thread [x−M,y−M,z−M] =

getNewValue((x,y,z),padded,Z)

ﬁnalOutput = np. zeros ( output .shape) #output image

for x in range(0, output .shape[0]) :

for y in range(0, output .shape[1]) :

for z in range(0, output .shape[2]) :

if not (x,y,z) in centers :

estimator = output[x−A:x+A+1,

y−A:y+A+1,

return output thread

def parallelize voxel selection (x):

output thread = np. zeros ( data .shape)

for y in range(M,padded.shape[1]−M):

for z in range(M,padded.shape[2]−M):

output thread [x−M,y−M,z−M] =

getNewValue voxel((x,y,z),padded)

return output thread

z−A:z+A+1]

count =

np.count nonzero( estimator )

sumOut = np.sum(estimator)

div = 0

if (count != 0):

div = sumOut/count

ﬁnalOutput [x,y,z] = div

else :

Passing to the functions shown before are just

values of X, this is because in this implementation

the chosen sub volume is a slice of the overall

input image, which is deﬁned as one value of

X and then within the processes the non local

means will be applied to the slice for the whole

of the Y and Z axis for that value of X. How

this is created is shown below with the call to

the classical implementation, the voxel selection

follows a similar pattern.

ﬁnalOutput [x,y,z] = 0

addition = np. zeros ( output .shape) #output image

addition = np.add(output , ﬁnalOutput )

return addition

D. Multi-Processor

Applying multi-threading to the aforementioned

algorithms allowed for full utilization of the com-

putational power available and created a drastic

start time classic thread = time . time ()

p classic = Pool(processes=3)





res classic =

p classic .map( parallelize classic

,range(M,padded.shape[0]−M))

p classic . close ()

p classic . join ()

print (”seconds ( classical ( thread )): ” ,

(time . time () − start time classic thread ))

for i in range( data .shape[0]) :

output classic thread [i]

\=

np.array ( res classic [i ][ i ])

Fig. 2. Noisy

Here there are three processes in the pool that

are available for the values of X that are passed

to it with in the .map() call this function call also

has a reference to the function deﬁnition that is

running the classic implementation and then after

the processes are done there are two calls one to

close() and one to join() in order to prevent mem-

ory leaks within Python. Now the return values

from each of the processes have to be unpacked

into the ﬁnal image which is what is occurring in

the loop with the values in resClassic at each layer

are loaded into the ﬁnal output image.

As we can see, the second image is much

more noisy than the ﬁrst one. Moreover, in the

next sections we will show how the results vary

by changing the size of the search volume. It

is important to notice that due to computational

limitation the ﬁlter is applied only to 3 layer of the

3D image therefore a degradation of the results is

expected for large values of M since more zeros

will be taken into consideration for the Euclidean

distance.

V. RESULTS

A. Classical

For our result section, we have processed the

original image taken form Brainweb database and

added Gaussian noise to it. The image that will

be processed with NL-mean algorithm is the one

with added Gaussian noise and the resulting image

will be compared using PSNR as a metric with the

original image. The difference between the two is

shown in ﬁgure below.

Fig. 3. Original-M1

Fig. 1. Ground





B. Classical with Voxel Selection

Fig. 4. Original-M2

Fig. 6. Voxel Selection-M1

Fig. 5. Original-M3

Fig. 7. Voxel Selection-M2

As expected the image is getting blurred with

increasing M. In the table below the resulting

PSNR results are shown.

M

1

PSNR (dB)

33.6

2

36.7

3

44.8

As proven by table above, the difference be-

tween ground truth and ﬁltered image increase

with an increase in search volume. This is coherent

with theoretical expectations.

Fig. 8. Voxel Selection-M3





The PSNR results are summarized in table be-

low.

M

1

PSNR (dB)

36.5

2

36.6

3

36.6

As it can be seen from the table above, there is

no degradation in the PSNR with an increase in M.

This is because the voxels that have 0 intensities

are disregarded by the selection process of the

most similar voxel. The main issues as it will

be shown in the multiprocessing section is that

of computational time. That is, due to a wrong

implementation of this improvement we did not

Fig. 11. BlockWise-M3

The PSNR results for Blockwise are shown in

see a decrease in computational time when voxel table below.

selection is used.

M

1

PSNR (dB)

48.3

C. Blockwise

2

65.3

3

62.3

As seen in the table, the PSNR results for

the Blockwise implementation increase as the M

increases and also deviate from the papers results

of around 30. There are a number of factors for

this, one is due to the abstraction nature of Block-

wise that is implemented in this report. With the

increase of the neighborhood size the M, the blocks

are then compared to blocks that are far away

from themselves, without voxel selection and fully

disregarding non similar blocks, regular Blockwise

will still factor these blocks in. Another factor that

can be attributed to the increase of the PSNR when

the M increases is the unpacking of the voxels

stage of blockwise. When M is increased and

the neighborhood becomes larger more smoothing

occurs and the PSNR increases.

Fig. 9. BlockWise-M1

D. Multiprocessing

Following are the graphs and tables, showing

the time complexity impact comparing the single

processor implementations to their multiprocessor

counterparts.

M

1

2

Classic time(s)

268.13

1248.77

Classic threaded time(s)

93.14

433.35

1182.2

Fig. 10. BlockWise-M2

3

3397.15





added to the time it took to run. In terms of re-

sults however, as mentioned before voxel selection

allowed for better edge preservation.

VI. PROBLEMS FACED

There were a few problems that were faced

while implemented the core features of this paper,

they ranged from issues with computation power

available, misunderstandings with how python ma-

trices index, or simply issues with incorrect initial

implementations.

Fig. 12. Classic Computation Time-M1

A. Size of input image

First with the classic implementation the impact

of the multiprocessing has a signiﬁcant impact on

the time to run the algorithm when M is increased.

When M = 3 we see a speedup of nearly 3

which makes a huge difference when running the

algorithm on a large image.

The ﬁrst hurdle that was faced in implementa-

tion of the algorithm was that the complexity of

the classic non local means makes it impossible

to make a slight change in the implementation,

run the algorithm on the whole image and see the

output. In order to get around this size complexity

issue, far smaller images were used in the ﬁrst

testing stages and often only a few layers were

selected for the algorithm to run on when wanting

to see an entire slice output.

M

1

Voxel time(s)

904.16

Voxel threaded time(s)

315.14

2

4218.70

1441.41

3

11425.09

3931.50

B. Selection of correct normalizing constant Z

While for the classical implementation of NL-

mean the normalizing constant can be know a

priori from the formula mentioned in section 4.A.

The same can not be done with the voxel selection

implementation. That is, the normalizing constant

can’t be known a priori because the number of vox-

els that will be selected for the calculation of the

Euclidean distance can not be known beforehand.

For that reason, the function ZVoxelSelection is be-

ing called in GetNewValueVoxel function across the

whole search volume so that the total Z is returned

and then passed to WVoxelSelection where it is

Fig. 13. Voxel Selection Computation Time -M1

With the multiprocessor improvement applied to used too normalize the exponential function. This

the voxel selection the speedup at M =3 was implementation created a huge computation over-

slightly lower than the speedup seen in classical, head that increased the computational time with

it is around 2.75. This however is still a impactful respect to classical implementation rather than

speedup considering the complexity of the algo- reducing it. A possible solution to this problem is

rithm. Comparing the times for the voxel selection to perform only 3 for loops in GetNewValueVoxel

verses the classic, voxel selection took longer due rather than 6. In fact, the position of voxel selection

to the fact that the voxel decision calculations could be returned directly by ZVoxelSelection so

happen dynamically while the program is running, that there will be no need to go across the whole

and thus adding to the complexity which in turn search volume one more time.





C. Python negative indexing

blockwise and the multi-threading implementa-

tions. The implementations of the classical al-

gorithm with smoothing parameters seems to be

consistent with theory expectations.The implemen-

tation with voxel selection is correct in the sense

that the PSNR does not vary with changing M

since the voxels with 0 intensity are excluded from

the selection. On the other hand, the computational

time for this improvements is large with respect

to the expected theoretical reduction. This is due

to a wrong implementation of the method for

ﬁnding the normalizing factor Z. The Blockwise

implementation caused many issues with overhead

and loss of complexity in the output, however in

the end with the ﬁnal implementation produced a

During the iterative process of building the code

implementation, it was noticed that python allows

for negative indexing of matrices. meaning that

passing a negative value into the matrix causes it

to return a value on the other side, to account for

the and to not allow the neighborhood aspect of

the non local means algorithm to get to a negative

index the solution was to pad the image to the size

of the search volume ensuring that even voxels of

the edges of the image wouldn’t cause a negative

indexing complication.

D. Blockwise

There were a few issues that occurred during good visual result with a time speedup compared

the implementation of the block wise improve- to classical. Applying the last improvement of

ment.One of which was an issue of overhead multiprocessing allowed for a drastic improvement

causing the algorithm to slow down. The ﬁrst upon the time complexity of the classical and voxel

attempt at the Blockwise implementation was to selection and is therefore considered a successful

use data structures to hold the locations of the implementation. The full code can be found in

centers of the blocks, all the voxels location’s [GitHub](https://github.com/EwanMcNeil/ELEC-444-PROJ.git)[ ](https://github.com/EwanMcNeil/ELEC-444-PROJ.git)and the main function to run is called

stored within them, and the block values that come main.py.

out of the non local means algorithm. This however

Some possible suggestions to tune the hyper-

caused a massive unnecessary overhead with too parameters of the algorithms such as M, b, s1,

many values beaning stored for no purpose. To ﬁx u1 and d is to use a stochastic optimizer were

this issue the block centers were known due to the objective function is to minimize the value of

having the N stepping size set and then inferring PSNR between ﬁltered image and the ground truth.

all the voxels overlapping within the block by

VIII. REFERENCES

using the volume deﬁned by the A value. Another

issue with Blockwise was that our implementation

had difﬁculty being used in the multi threading

since the implementation of Blockwise followed a

more sequential steps that relied upon each other, it

became difﬁcult to implement the multi-threading.

Comparatively however the time complexity of a

classical implementation verses Blockwise imple-

mentation still has a large speedup.

[1] ”An optimized blockwise nonlocal

means denoising ﬁlter for 3-D magnetic

resonance images” IEEE transactions on

medical imaging vol. 27,4 (2008): 425-41.

doi:10.1109/TMI.2007.906087

VII. CONCLUSION

In conclusion, in this paper we have imple-

mented the NL-means algorithm and all the im-

provements proposed by Coupe et al. (2008) to

deal with the huge computational complexity in-

herently present in the classical algorithm. That

is we have implemented the automatic toning

of the smoothing parameter, the voxel selection

to pre-select only the most similar voxels, the


