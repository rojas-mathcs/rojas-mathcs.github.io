# Copyright 2013, J. Maurice Rojas
# Modified by Matthew J. Barry for Homework 4

# utility function to split a list {l} into chunks of size {size}
def chunks(l, size):
    for i in xrange(0, len(l), size):
        yield l[i:i+size]

# got these primes from http://primes.utm.edu/lists/small/
# both have about 130 digits
p = 49105133753678962419391374496727457097027115514864975959094375339860082156189031654878635228598019107441653554142227931036658923353143622433
q = 64245974067698974727308576291067189905635129662015537917175582385219809100094736540969262037560193877488366836626125921856719528308413429711

# could've picked smaller p and q: all we need is n with >138 digits
n = p*q

# encryption / decryption keys
e = 9872398472439
d = 1/e % (p-1)*(q-1)

# load full plaintext of Raymond T. Chandler's **The Big Sleep** as list of numbers: 'a'=>0, 'b'=>1, etc.
nums = map(lambda c: ord(c)-ord('a'), open('plaintext.txt', 'r').read().strip())

# split nums into {csize}-sized chunks
csize = 97

# treat elements of chunk as digits of a base-26 number
base = [26^i for i in xrange(csize)] 
newnums = (reduce(lambda x,y: x+y, map(lambda x,y: x*y, base, c)) for c in chunks(nums, csize))

# raise chunk-number to {e}-th power (mod {n}) to encrypt (offset by 2 since 0^k=0 and 1^k=1 for all k)
ciphernums = (Mod(num+2, n)^e for num in newnums)

# convert back to list of characters once encrypted
ciphertext = [chr(c+ord('a')) for cnum in ciphernums for c in Integer(cnum-2).digits(26)]

# Perform frequency analysis on ciphertext
freqs = {}
for char in ciphertext:
    if char in freqs:
        freqs[char] += 1
    else:
        freqs[char] = 1

# output results
for l in sorted([(-v,k) for k,v in freqs.iteritems()]):
    print("%c: %5.3f" % (l[1], -float(l[0])/len(ciphertext)))

