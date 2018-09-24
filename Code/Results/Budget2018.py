import numpy as np

dollarRandConv = 13.5
myTtlBudget = 22807.*dollarRandConv
myRate = 300.
#  3*4*7*2  #the no hours I should work in a quarter
my2017Invoices = np.array([24300, 41700, 3*4*7*2.*myRate])   # last amount is a guess for Dec
# fixedCosts2018 = 30*32.50*14.5      # 30km2 x $32.50/km2   # worst case costs for imagery
fixedCosts2018 = 0.                 # we got a grant from dg, so we can do something with the time now
varCosts2018 = myTtlBudget - my2017Invoices.sum() - fixedCosts2018
hoursPerDay = 7.

print "Total budget: ", myTtlBudget
print "2017 invoices: ", my2017Invoices.sum()
print "2018 remainder: ", myTtlBudget - my2017Invoices.sum()
# print "2018 remainder in days: ", (myTtlBudget - my2017Invoices.sum())/(hoursPerDay*myRate)
# print "2018 remainder in weeks: ", (myTtlBudget - my2017Invoices.sum())/(hoursPerDay*myRate*5)
print "2018 fixed (imagery) costs: ", fixedCosts2018
print "2018 var (time worked) costs: ", varCosts2018
print "2018 time: %f (hrs), %f (days), %f (wks), %f (months)" % (varCosts2018/myRate, varCosts2018/(myRate*hoursPerDay), varCosts2018/(myRate*hoursPerDay*2), varCosts2018/(myRate*hoursPerDay*2*4))

# how to split time over quarters?
# daysPerWeekForQtr = np.array([1., 2., 2.5, 2.5])
# daysPerWeekForQtr = np.array([1., 1., 1.5, 4.5])
daysPerWeekForQtr = np.array([1., 1., 1.5, 1.5, 1.5, 2.])        # extend to 2019
(daysPerWeekForQtr*hoursPerDay*4*3).cumsum()
(daysPerWeekForQtr*hoursPerDay*4*3).cumsum() - varCosts2018/myRate   #in which quarter do we complete the ttl hours?

# how many months do we need with the above time split
# (daysPerWeekForQtr*7*4*3).cumsum() - varCosts2018/myRate   #in which quarter do we complete the ttl hours?
idx = np.where(((daysPerWeekForQtr*hoursPerDay*4*3).cumsum() - varCosts2018/myRate) > 0)[0]   #in which quarter do we complete the ttl hours?

tmpHours = (varCosts2018 / myRate) - (daysPerWeekForQtr*hoursPerDay*4*3).cumsum()[idx-1]   #hours to be worked in final quarter
tmpDays = tmpHours/(hoursPerDay)   #hours to be worked in final quarter
tmpWeeks = tmpHours/(daysPerWeekForQtr[idx]*hoursPerDay)   #weeks to be worked in final quarter
ttl = 0
for i, d in enumerate(daysPerWeekForQtr):
    if i<idx:
        print "Days to be worked in qtr %i: %f, 5 day weeks: %f" % (i + 1, d*4.*3., d*4.*3./7.)
    elif i==idx:
        print "Days to be worked in qtr %i: %f, 5 day weeks: %f" % (i + 1, tmpDays, tmpDays/7.)
    else:
        break

print "Weeks to be worked in qtr %i: %f" % (idx+1, tmpWeeks)

budgetPerQuarter = daysPerWeekForQtr*hoursPerDay*myRate*4*3
budgetPerQuarter[-1] = tmpHours*myRate
print 'var costs budgetPerQuarter: ', budgetPerQuarter
print '$ var costs budgetPerQuarter: ', budgetPerQuarter/dollarRandConv
budgetPerQuarter[-1] += fixedCosts2018   # buy imagery in first quarter but only invoice in last quarter
print 'all costs budgetPerQuarter: ', budgetPerQuarter
print '$ all costs budgetPerQuarter: ', budgetPerQuarter/dollarRandConv

budgetPerQuarter.sum()  #check we sum to "2018 remainder"

# how many days can we fund with the remaining budget?
(budgetPerQuarter.sum() - budgetPerQuarter[:3].sum())/(hoursPerDay*myRate)