# SPDZ-2-With Extensions for Ring Supporting a Floating-point Operation 

A fork of the University Of Bristol SPDZ-2 Repository, with changes to support extending the SPDZ-2 Framework to run additional protocols. Changes performed by Bar Ilan Cryptography Research Group and NEC Security Research Labs.
This repository is the improved one of SPDZ-2(https://github.com/nec-mpc/SPDZ-2). In this version, the secure computation on floating point numbers via ring-based 3PC protocol is newly supported.

We would like to thank to the team behind the SPDZ-2 framework, which is an extensive effort and an excellent contribution to the MPC community. Special thanks to Marcel Keller for his numerous insights and explanations making this work possible.

(C) 2017 University of Bristol. See License.txt Software for the SPDZ and MASCOT secure multi-party computation protocols. See Programs/Source/ for some example MPC programs.

## SPDZ-2 With Extensions - rationale

The SPDZ-2 extensions is a mechanism that enables substitution of the original implementation of various operations with an alternate external implementation. This is done by dynamically loading a configured library and prescribed API function pointers. In runtime, the SPDZ-2 processor will call the loaded API functions instead of the original implementation and provide it with the required parameters. In this repository, we set library for ring based protocol as configured library.

### MPC programs source code
The [Programs/Source](https://github.com/nec-mpc/SPDZ-2-Float/tree/master/Programs/Source) folder of this fork contains MPC programs added as part of our work to evaluate different protocols under the framework. For example, the following program performs some basic operations.  
```
##### examples for secret integer operation #####

# sint: secret integer
# cint: clear integer

## addition ##

a1 = sint(1)
b1 = sint(2)
c1 = a1 + b1
# expected value: c1=3
print_ln("c1=%s", c1.reveal())

a2 = sint(2)
b2 = cint(1)
c2 = a2 + b2
# expected value: c2=3
print_ln("c2=%s", c2.reveal())


## multiplication ##

# "sint * sint" requires communications. 
a3 = sint(1)
b3 = sint(2)
c3 = a3 * b3
# expected value: c3=2
print_ln("c3=%s", c3.reveal())


# "sint * cint" does not require communications.
# It needs only local operations. Hence, it takes no communication cost.
a4 = sint(2)
b4 = cint(1)
c4 = a4 * b4
# expected value: c4=2
print_ln("c4=%s", c4.reveal())

## comparison ##
a5 = sint(5)
b5 = sint(3)
c5 = (a5 < b5)
# expected value: c5=0
print_ln("c5=%s", c5.reveal())

a6 = sint(4)
b6 = sint(10)
c6 = (a6 < b6)
# expected value: c6=1
print_ln("c6=%s", c6.reveal())

## equality testing ##
a7 = sint(1)
b7 = sint(2)
c7 = (a7 == b7)
# expected value: c7=0
print_ln("c7=%s", c7.reveal())

a8 = sint(5)
b8 = sint(5)
c8 = (a8 == b8)
# expected value: c8=1
print_ln("c8=%s", c8.reveal())

## secure if_else ##
# Secure if_else operation is branching program while hiding the condition value.
# Let cond be the condition value which is 1 or 0.
# For example, cond.if_else(x, y) outputs x, otherwise it ouputs y.
# That is, cond.if_else(x, y) equals to cond * x + (1 - cond) * y.

a9 = sint(10)
b9 = sint(8)
cond1 = (a9 < b9)
c9 = cond1.if_else(100, 20)
# expected value: c9=20
print_ln("c9=%s", c9.reveal())

a10 = sint(10)
b10 = sint(10)
cond2 = (a10 == b10)
c10 = cond2.if_else(100, 20)
# expected value: c10=100
print_ln("c10=%s", c10.reveal())

##### examples for secret fixed operation #####

## sfix: secure fixed value
## cfix: clear fixed value 

## addition ##
a11 = sfix(1.5)
b11 = sfix(3.21)
c11 = a11 + b11
# expected value: c11=4.71
# actual output: c11=4.710007
print_ln("c11=%s", c11.reveal())

a12 = sfix(1.5)
b12 = cfix(0.21)
c12 = a12 + b12
# expected value: c12=1.71
# actual output: c12=1.710007
print_ln("c12=%s", c12.reveal())

## multiplication ##
a13 = sfix(1.5)
b13 = sfix(0.2)
c13 = a13 * b13
# expected value: c13=0.3
# actual output: c13=0.300003
print_ln("c13=%s", c13.reveal())

a14 = sfix(1.7)
b14 = cfix(0.3)
c14 = a14 * b14
# expected value: c14=0.51
# actual output: c14=0.510010
print_ln("c14=%s", c14.reveal())

## division ##
a15 = sfix(3.5)
b15 = sfix(2.1)
c15 = a15 / b15
# expected value: c15=1.6666...
# actual output: c15=1.666656
print_ln("c15=%s", c15.reveal())


##### examples for secret float operation #####

## sfloat: secure float value

## addition ## 
a16 = sfloat(3.2)
b16 = sfloat(1000000000.1)
c16 = a16 + b16
# expected value: c16=1000000003.3
# actual output: c16=1.000000e+09
print_ln("c16=%s", c16.reveal())
 
## multiplication ##
a17 = sfloat(0.00000001)
b17 = sfloat(1000000000)
c17 = a17 * b17
# expected value: c17=10
# actual output: c17=9.999999e+00
print_ln("c17=%s", c17.reveal())

## division ## 
a18 = sfloat(1.2345678)
b18 = sfloat(0.0001)
c18 = a18 / b18
# expected value: c18=12345.678
# actual output: c18=1.234568e+04
print_ln("c18=%s", c18.reveal())

print_ln("LT Test:")
x = sfloat(2*(10**5))
y = sfloat(2*(10**(-5)))
res = (x<y)
print_ln("res == %s", res.reveal())
```
### SPDZ-2 extension library
See https://github.com/nec-mpc/SPDZ-2-Extension-Ring-Float for our implemented extension library.

## SPDZ-2
### Requirements:

- GCC (tested with 4.8.5) or ICC (18.0.3)
- MPIR library, compiled with C++ support (use flag --enable-cxx when running configure)
- libsodium library, tested against 1.0.11
- CPU supporting AES-NI
- Python 2.x (tested with 2.7.5)

### To compile:
1) Download files of this repository to above environment.

2) Change directories to download one.

3) Run `make clean all`

### To generate the bytecode:
1) Set the program source file of MPC on [Programs/Source](https://github.com/nec-mpc/SPDZ-2-Float/tree/master/Programs/Source). 
 - File extension is ".mpc"
 
2) Change directories to download one.

3) Run `python compile.py [PROGRAM NAME]`

### To run the protocol:
1) Set environment variables for extension library. 
 - See the URL and README.md

2) Change directories to download one.

3) Run each entity as follows.
* [Proxy]
	`./Server.x 3 [port number]`
	
* [Each MPC server(0/1/2)] 
	`.Player-Online.x -pn [port number] -lgp 64 [server ID] [PROGRAM NAME]`

- test_inner_product (This program computes the inner-product value)
   * [Proxy]
     `./Server.x 3 60000`
   * [MPC server0]
     `.Player-Online.x -pn 60000 -lgp 64 0 test_inner_product`
   * [MPC server1]
      ` .Player-Online.x -pn 60000 -lgp 64 1 test_inner_product`
   * [MPC server2]
      `.Player-Online.x -pn 60000 -lgp 64 2 test_inner_product`
