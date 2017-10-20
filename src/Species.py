import re

## species map
species = {}

class Species:
  def __init__(self, spString):
  
    ## parse the species for charge and state
    def parseSpecies(str):
      m = re.match('^\((.*)\)(?:(([1-9]?|[1-9][0-9]*)([\+-]))|([lgs]|aq))$', str)
      
      if m is None:
        raise SystemExit("Species %s is ill-formed (state)" % str)
      
      sp = m.group(1)
      
      ## if there is a charge, get charge
      if m.group(2) is not None: 
        (mag, sign) = (m.group(3), m.group(4))
      
        if mag == "":
          charge = 1
        else:
          charge = int(mag)
      
        if sign == "-":
          charge *= -1
        
        ## charged species must be in solution        
        state = "aq"
    
      ## if no charge, then get state
      else:
        charge = 0
        state = m.group(5)
        
      
      ## parse the formula to get element totals
      def parseComposition(component, n):
        if component == "":
          return {}
      
        ## flag for parenthesized component
        pFlag = False
      
        ## if next component is in brackets
        ##   we must manually find what is between brackets
        ##   (regex can't do this as it has no state)
        if component[0] == '(':
          parens = 1
          l = 0
          
          for c in component[1:]:
            if c == '(':
              parens += 1
            elif c == ')':
              parens -= 1
            
            if parens == 0:
              break
            
            l += 1
          
          ## report error if unbalanced parens, or nothing between parens
          if parens != 0 or l == 0:
            raise SystemExit("Species %s is ill-formed (parentheses)" % str)
          
          ## take what is between parens          
          c = component[1:1 + l] 
          
          ## regex to find number of occurrences and rest of formula
          m = re.match('^([1-9][0-9]*)?(.*)$', component[2 + l:]) 
          (num, rest) = (m.group(1), m.group(2))
          
          pFlag = True
        
        ## otherwise just use full regex
        else:                
          m = re.match('^([A-Z][a-z]?)([1-9][0-9]*)?(.*)$', component)
        
          if m is None:
            raise SystemExit("Species %s is ill-formed (formula)" % str)
            
          (c, num, rest) = (m.group(1), m.group(2), m.group(3))
        
        
        dictRest = {}
        
        ## if there is more to the formula, 
        ##   then parse it and store totals in dictRest
        if rest is not None:
          dictRest = parseComposition(rest, n)
        
        
        if num is not None:
          dictComp = parseComposition(c, n * int(num))
        ## if there is no number and component was parenthesized
        ##   then there is still parsing to do
        elif pFlag: 
          dictComp = parseComposition(c, n)
        else:
          dictComp = { c: n }
        
        ## add the total of the current component 
        ##   into the dictionary  
        for key, value in dictComp.items():
          dictRest[key] = dictRest.get(key, 0) + value
        
        return dictRest
        
      comp = parseComposition(sp, 1)
          
      return (state, charge, comp)
    
    ## if we've stored this species already, just retrieve it
    if spString in species:
      sp = species[spString]
      self.formula = sp.formula
      self.state = sp.state
      self.charge = sp.charge
      self.composition = sp.composition
      return
    
    ## otherwise, parse and store in species    
    self.formula = spString   
    (self.state, self.charge, self.composition) = parseSpecies(spString)
    species[spString] = self