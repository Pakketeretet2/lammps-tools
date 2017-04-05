import dumpreader
import sys
import time

dump_name = "sample_strip.dump"
d = dumpreader.dumpreader( dump_name )
b = d.getblock()

while not d.at_eof:
    print("t = ", b.meta.t)
    b = d.getblock()

d  = dumpreader.dumpreader( dump_name )
b  = d.getblock()
dumpreader.block_to_dump_write( b, "/dev/stdout", "w" )

d2 = dumpreader.dumpreader_cpp( dump_name )
b = d2.getblock()

dumpreader.block_to_dump_write( b, "/dev/stdout", "w" )


print("Testing skip speed.",file=sys.stderr)
d = dumpreader.dumpreader( dump_name )
skip = 55


d = dumpreader.dumpreader( dump_name )
t1 = time.clock()
for i in range(0,skip):
    b = d.getblock()
b = d.getblock()
t2 = time.clock()
print("Dumb python:  ", t2-t1, " s, tstep = ", b.meta.t)


d = dumpreader.dumpreader( dump_name )
t1 = time.clock()
d.fast_forward(skip)
b = d.getblock()
t2 = time.clock()
print("Smart python: ", t2-t1, " s, tstep = ", b.meta.t)



d = dumpreader.dumpreader_cpp( dump_name )
t1 = time.clock()
for i in range(0,skip):
    b = d.getblock()
b = d.getblock()
t2 = time.clock()
print("Dumb C++:     ", t2-t1, " s, tstep = ", b.meta.t)

d = dumpreader.dumpreader_cpp( dump_name )

t1 = time.clock()
d.fast_forward(skip)
b = d.getblock()
t2 = time.clock()
print("Smart C++:    ", t2-t1, " s, tstep = ", b.meta.t)

