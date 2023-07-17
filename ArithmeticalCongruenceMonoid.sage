# Class: ArithmeticalCongruenceMonoid
# 
# Authors: Jacob Hartzer and Christopher O'Neill
# Version: 0.1a
# 
# Last Edited 5/18/2016
# 
# This is a python class based around the factorization of Arithmetical Congruence 
#   Monoids (ACMs). Most of the operations are primitive, but nonetheless useful
#   for investigating the factorization structure of ACMs. 
# 
# Feel free to copy and make changes, but please leave the original (Dropbox) version unedited.  
# More functionality will be added, and feature requests are welcome, just email the authors
# 
#####################################################################################
# 
# Usage
# 
# load('/path/to/file/ArithmeticalCongruenceMonoid.sage)
# 
# hilbert = ArithmeticalCongruenceMonoid(1,4)
# hilbert.Factorizations(441)
# hilbert.DeltaSet(441)
# 
#####################################################################################

import itertools
import collections import Counter

def IntegerDivisors(num):
	# Builds lists of possible divisors by all possible combinations of each prime factor
	factorization = factor(num)
	ret = [1]
	
	for i in range(len(factorization)):
		B = []
		for j in range(1,factorization[i][1]+1):
			for k in ret:
				B.append(k*(factorization[i][0])^j)
		ret = ret + B
	
	return ret

def DeltaListFromList(l):
	return [l[i+1] - l[i] for i in range(len(l)-1)]



class ArithmeticalCongruenceMonoid:
	
	def __init__(self, a, b):
		if ((a^2) % b) != (a % b):
			raise ValueError("a and a^2 must be equivalent mod b")
		
		self.b = b
		self.a = a % b
		
		self.__factorizations = {1:[[]]}
		self.__irreducibles = {1:False}
		self.__lengthsets = {1:[0]}

	def __closedDivisors(self, num):
		return [i for i in IntegerDivisors(num) if i in self and num/i in self]
	def Factorizations(self, num):
		if num in self.__factorizations:
			return self.__factorizations[num]
		divisors = sorted(self.__closedDivisors(num))
		for d in divisors:
			if d in self.__factorizations:
				continue
			self.__factorizations[d] = []
			for s in divisors:
				if(s >= d):
					break
				if(self.__irreducibles[s] and d/s in self):
					for f in self.__factorizations[(d/s)]:
						if(s >= f[-1]):
							self.__factorizations[d].append(f + [s])
			if(len(self.__factorizations[d]) == 0):
				self.__factorizations[d] = [[d]]
				self.__irreducibles[d] = True
			else:
				self.__irreducibles[d] = False
		return self.__factorizations[num]
	def NumberOfFactorizations(self, num):
		return len(self.Factorizations(num))

	def FactorizationsUpToElement(self, nmax):
		# Finds all arithmetic factorizations up to a certain element. Useful to run before long calculations
		for i in range(self.a, nmax + 1, self.b):
			self.Factorizations(i)
	
	def IsIrreducible(self,num):
		# checks if an element is irreducible. 
		if num in self.__irreducibles:
			return self.__irreducibles[num]
		
		if num == 1:
			return False
		
		if num in self.__factorizations:
			return max(self.LengthSet(num)) == 1
		
		# removes all numbers not divisors in the monoid
		divisors = [i for i in IntegerDivisors(num) if i in self and num/i in self]
		
		self.__irreducibles[num] = (len(divisors) <= 2)
		return self.__irreducibles[num]
	
	def IrreduciblesUpToElement(self,nmax):
		# creates a dict of the irreducibles up to a given element 
		start = self.a
		if self.a == 1:
			start = start + self.b
		
		# pre allocates dictionary memory
		for i in range(start,nmax,self.b):
			self.__irreducibles[i] = True
		
		# essentially a Sieve of Eratosthenes for ACMs
		for i in range(start,nmax,self.b):
			# for each product of two elements, sets dict(product) = 0
			for j in range(i,nmax/i+1,self.b):
				self.__irreducibles[i*j] = False
	
	def PlotNumberOfFactorizationsToElement(self,nmax):
		# plots number of factorizations up to a given element
		plot_list=[]
		for i in range(self.a+self.b,nmax+1,self.b):
			plot_list.append([i,self.NumberOfFactorizations(i)])
		return list_plot(plot_list)
	
	def LengthSet(self, num):
		if num in self.__lengthsets:
			return self.__lengthsets[num]
		divisors = sorted(self.__closedDivisors(num))
		for d in divisors:
			if d in self.__lengthsets:
				continue
			self.__lengthsets[d] = []
			for s in divisors:
				if( s >= d):
					break
				if(self.__irreducibles[s] and d/s in self):
					for l in self.__lengthsets[(d/s)]:
						self.__lengthsets[d].append(l + 1)
			if(len(self.__lengthsets[d]) == 0):
				self.__lengthsets[d] = [1]
				self.__irreducibles[d] = True
			else:
				self.__irreducibles[d] = False
			self.__lengthsets[d] = list(set(self.__lengthsets[d]))
		return self.__lengthsets[num]
	
	def DeltaSet(self, num):
		return list(sorted(set(DeltaListFromList(self.LengthSet(num)))))
	def LengthSet(self,num,t):
		if t == 1:
			return self.LengthSet(num)
		if t == 0:
			return set(len(set(F)) for F in self.Factorizations(num))
		if t == oo:
			return set(Counter(F).most_common()[0][1] for F in self.Factorizations(num))
		print('ERR: t value not implemented, defaulting to t = 1')
	def MaxFactorizationLength(self,num):
		# finds the max factorization length of a given element
		return max(self.LengthSet(num))
	
	def MinFactorizationLength(self,num):
		# finds the min factorization length of a given element
		return min(self.LengthSet(num))
	def MaxFactorizationLength(self,num,t):
		return max(self.LengthSet(num,t))
	def MinFactorizationLength(self,num,t):
		return min(self.LengthSet(num,t))
	def Elasticity(self,num):
		# finds elasticity of an element: max factorization length / min...
		if num == 1:
			return 1
		return float(self.MaxFactorizationLength(num))/float(self.MinFactorizationLength(num))
	
	@staticmethod
	def MonoidsForGivenB(b):
		return [a for a in range(b) if (a^2 % b) == a]
	
	def __eq__(self, other):
		return (self.a, self.b) == (other.a, other.b)
	
	def __ne__(self, other):
		return (self.a, self.b) != (other.a, other.b)
	
	def __contains__(self, other):
		return int(other) == other and (other == 1 or other % self.b == self.a)
	
	def __repr__(self):
		return "Arithmetical Congruence Monoid " + str((self.a, self.b))
	
	def __str__(self):
		return self.SaveToString()
	
	
	
	
	
	
	# def LongestFactorizationsInRange(self, nmax):
	# 	# returns indices and length of element with longest factorization
	# 	self.FactorizationsUpToElement(nmax)
	# 	index = []
	# 	max_length = 0
	# 	for i in range(self.a,nmax+1,self.b):
	# 		if len(self.Factorizations(i)) > max_length:
	# 			max_length = len(self.Factorizations(i))
	# 			index = [i]
	# 		elif len(self.Factorizations(i)) == max_length:
	# 			index.append(i)
	# 	return index, max_length
	
	# def NumberOfFactorizations(self,num):
	# 	# returns the number of factorizations of an element
	# 	return len(self.Factorizations(num))
	
	# def PlotMaxFactorizationLengthToElement(self,nmax):
	# 	# plots max factorization length up to a given element
	# 	plot_list=[]
	# 	for i in range(self.a+self.b,nmax+1,self.b):
	# 		plot_list.append([i,self.MaxFactorizationLength(i)])
	# 	return list_plot(plot_list)
	
	# def PlotMinFactorizationLengthToElement(self,nmax):
	# 	# plots min factorization length up to a given element
	# 	plot_list=[]
	# 	for i in range(self.a+self.b,nmax+1,self.b):
	# 		plot_list.append([i,self.MinFactorizationLength(i)])
	# 	return list_plot(plot_list)
	
	# def MaxElasticityToElement(self,nmax):
	# 	# Finds max elasticity up to given element.
	# 	max_elasticity = 0
	# 	for i in range(self.a+self.b,nmax+1,self.b):
	# 		if self.Elasticity(i) >= max_elasticity:
	# 			max_elasticity = self.Elasticity(i)
	# 	return max_elasticity
	
	# def PlotElasticitiesToElement(self,nmax):
	# 	# plots elasticities to a given element 
	# 	plot_list=[]
	# 	for i in range(self.a+self.b,nmax+1,self.b):
	# 		plot_list.append([i,self.Elasticity(i)])
	# 	p = list_plot(plot_list,figsize=(10,4),size = 5,title = 'Elasticities of M_4,6')
	# 	return p
	
	# def PlotIrreduciblesToElement(self,nmax):
	# 	# plots irreducibles to a given element (irreducible = 1)
	# 	self.IrreduciblesUpToElement(nmax)
	# 	p = list_plot(self.__irreducibles,figsize=(10,4),size = 50,title = 'Irreducibles of M_4,6',ymin = -0.5,ymax = 1.5)
	# 	return p
	
	# def FindMaxDeltaSet_Irreducibles(self,nmax):
	# 	# Finds the largest distance between two irreducible elements 
	# 	self.IrreduciblesUpToElement(nmax)
	# 	max_delta = 1
	# 	start = self.a
	# 	if start == 1:
	# 		start = start + self.b
	# 	hold = start
	# 	for i in range(start,nmax,self.b):
	# 		if self.__irreducibles[i] == 1:
	# 			if (i-hold)>max_delta:
	# 				max_delta = i-hold
	# 				print max_delta,i,factor(i)
	# 			hold = i
	
	# def FindMaxDeltaSet_Reducibles(self,nmax):
	# 	# Finds the largest distance between two reducible elements 
	# 	self.IrreduciblesUpToElement(nmax)
	# 	max_delta = 1
	# 	start = self.a
	# 	if start == 1:
	# 		start = start + self.b
	# 	hold = start 
	# 	for i in range(start,nmax,self.b):
	# 		if self.__irreducibles[i] == 0:
	# 			if (i-hold)>max_delta:
	# 				max_delta = i-hold
	# 				print max_delta,i,factor(i)
	# 			hold = i
	
	# def PlotDeltaSet_Reducibles(self,nmax):
	# 	# plots the delta set of reducible elements to a given element.
	# 	# delta set is the distance between consecutive reducibles
	# 	self.IrreduciblesUpToElement(nmax)
	# 	hold = self.a
	# 	delta_list = []
	# 	for i in range(self.a+self.b,nmax,self.b):
	# 		if self.__irreducibles[i] == 0:
	# 			delta_list = delta_list + [i-hold]
	# 			hold = i
	# 	return list_plot(delta_list)
	
	# def PlotDeltaSet_Irreducibles(self,nmax):
	# 	# plots the delta set of irreducible elements to a given element.
	# 	# delta set is the distance between consecutive irreducibles
	# 	self.IrreduciblesUpToElement(nmax)
	# 	hold = self.a
	# 	delta_list = []
	# 	for i in range(self.a+self.b,nmax,self.b):
	# 		if self.__irreducibles[i] == 1:
	# 			delta_list = delta_list + [i-hold]
	# 			hold = i
	# 	return list_plot(delta_list)
	
	# def PlotVariableDeltaSet_Reducibles(self,nmax,start,period):
	# 	# if you want to change the period in which the program looks for irreducibles
	# 	# this will do that for you. 
	# 	self.IrreduciblesUpToElement(nmax)
	# 	hold = self.a + start*self.b
	# 	delta_list = []
	# 	for i in range(self.a + start*self.b, nmax, period*self.b):
	# 		if self.__irreducibles[i] == 0:
	# 			delta_list = delta_list + [(i-hold)/period]
	# 			hold = i
	# 	return list_plot(delta_list)
	
	# def Percent_Irreducible(self,nmax):
	# 	# returns the percent of elements that irreducible to a given element 
	# 	self.IrreduciblesUpToElement(nmax)
	# 	return float(sum(self.__irreducibles.values())/len(self.__irreducibles))
		
