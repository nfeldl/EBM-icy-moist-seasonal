#!/software/anaconda3/bin

cday    = 86400.0         # sec in calendar day ~ sec
cyear   = cday*365.25     # sec in calendar year ~ sec
cpsw    = 3.996e3         # specific heat of sea water ~ J/kg/K
rhosw   = 1.026e3         # density of sea water ~ kg/m^3

def heatcapacity(h): 

  """
  Calculate mixed-layer heat capcity in units 
  W yr m^-2 K^-1 given depth in m
  """

  heat_capacity = cpsw*rhosw*h 

  heat_capacity = heat_capacity / cyear

  return heat_capacity

def depth(C): 

  """
  Calculate mixed-layer depth given heat capcity
  in units of W yr m^-2 K^-1
  """
  C = C * cyear # convert to Joules

  depth = C/(cpsw*rhosw)

  return depth 

def main():

  print(f'Mixed layer depth is {depth(9.8)}')
 
  print(f'Mixed layer heat capacity is {heatcapacity(30)}')

if __name__ == '__main__':
  main()



