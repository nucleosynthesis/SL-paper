import model 

nb = len(model.data)

f1 = open("Table_1.yaml","w")

# Do a table for signal and data
f1.write("dependent_variables:\n")
f1.write("- header:\n")
f1.write("    name: Data Events\n")
f1.write("  values:\n")

for x in model.data: 
  f1.write("    - value: %d\n"%x)

f1.write("- header:\n")
f1.write("    name: Signal Events\n")
f1.write("  values:\n")

for x in model.signal: 
  f1.write("    - value: %g\n"%x)

f1.write("independent_variables:\n")
f1.write("- header:\n")
f1.write("    name: Bin Number\n")
f1.write("  values:\n")

for x in range(nb): 
  f1.write("  - value: %d\n"%(x+1))

f1.close()

# now do the table for m_1

f1 = open("Table_2.yaml","w")

f1.write("dependent_variables:\n")
f1.write("- header:\n")
f1.write("    name: $m_1$\n")
f1.write("  values:\n")

for x in model.background: 
  f1.write("  - value: %g\n"%x)

f1.write("independent_variables:\n")
f1.write("- header:\n")
f1.write("    name: Bin Number\n")
f1.write("  values:\n")

for x in range(nb): 
  f1.write("  - value: %d\n"%(x+1))

f1.close()

# and a table for m2 

f1 = open("Table_3.yaml","w")

f1.write("dependent_variables:\n")

ycounter = 0 
for x in range(nb): 
 f1.write("- header:\n")
 f1.write("    name: Bin %d\n"%(x+1))
 f1.write("  values:\n")
 for y in range(nb): 
   f1.write("  - value: %g\n"%(model.covariance[ycounter]))
   ycounter+=1
 
f1.write("independent_variables:\n")
f1.write("- header:\n")
f1.write("    name: Bin Number\n")
f1.write("  values:\n")
for x in range(nb): 
  f1.write("  - value: %d\n"%(x+1))

f1.close()



# and a table for m3

f1 = open("Table_4.yaml","w")

f1.write("dependent_variables:\n")
f1.write("- header:\n")
f1.write("    name: $m_3$\n")
f1.write("  values:\n")

for x in model.third_moment: 
  f1.write("  - value: %g\n"%x)

f1.write("independent_variables:\n")
f1.write("- header:\n")
f1.write("    name: Bin Number\n")
f1.write("  values:\n")

for x in range(nb): 
  f1.write("  - value: %d\n"%(x+1))

f1.close()

