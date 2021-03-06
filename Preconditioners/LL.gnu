a(q) = (2*(q/qf) + (q/qf)**2)**2 + d**2
b(q) = (2*(q/qf) - (q/qf)**2)**2 + d**2

c(q) = (2*(q/qf) - (q/qf)**2)/d
e(q) = (2*(q/qf) + (q/qf)**2)/d

tanterm(q) = atan( a(q) ) + atan( b(q) )
logterm(q) = log( e(q) ) - log( c(q) )

qf = 1
d=1


llq(q) = 1 + (2/(qf*3.14159)) * (   1/((q/qf)**2) - (d/(2*(q/qf)**3)) * (atan( c(q) ) + atan( e(q) )) +    (d**2/(8*(q/qf)**5) + 1/(2*(q/qf)**3) - 1/(8*(q/qf))) * (log( a(q) ) - log( b(q) ))  )


aa(q) = (2*(q/qf) + (q/qf)**2)**2 + dd**2
bb(q) = (2*(q/qf) - (q/qf)**2)**2 + dd**2

cc(q) = (2*(q/qf) - (q/qf)**2)/dd
ee(q) = (2*(q/qf) + (q/qf)**2)/dd

tanterm1(q) = atan( aa(q) ) + atan( bb(q) )
logterm1(q) = log( ee(q) ) - log( cc(q) )


dd=0.4

llq2(q) = 1 + (2/(qf*3.14159)) * (   1/((q/qf)**2) - (dd/(2*(q/qf)**3)) * (atan( cc(q) ) + atan( ee(q) )) +    (dd**2/(8*(q/qf)**5) + 1/(2*(q/qf)**3) - 1/(8*(q/qf))) * (log( aa(q) ) - log( bb(q) ))  )






set term pdf
set output "LL.pdf"

set xrange [0.2:4]


kerA = 0.6
kerB = 1.5

kerkerq(q) = (kerA * ((q)**2))/(kerB**2 + (q)**2)




limit(q) = 1 + (2/(qf*3.14159)) * ( 1/((q/qf)**2) )

#plot (1/llq(x)), 1/kerkerq(x), 1/limit(x)

#p 1/kerkerq(x)
p (1/llq(x)), (1/llq2(x)), kerkerq(x)
