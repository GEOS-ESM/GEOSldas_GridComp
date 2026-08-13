import numpy as np  

def is_leap_year(year=None):
 
    if year%400 == 0 or (year%4 ==0 and year%100 != 0):
        return True
    else:
        return False
    
def pentad_of_year(day_of_year=None,year=None):
    
    if is_leap_year(year) and day_of_year >= 59:
        pentad=np.floor((day_of_year - 2) / 5) + 1
    else:
        pentad=np.floor((day_of_year - 1) / 5) + 1
    
    return pentad

def days_in_month(year=None,month=None):

    days_in_month_leap=[31,29,31,30,31,30,31,31,30,31,30,31]
    days_in_month_nonleap=[31,28,31,30,31,30,31,31,30,31,30,31]
    if is_leap_year(year):
        n_days=days_in_month_leap(month)
    else:
        n_days=days_in_month_nonleap(month)

def get_dofyr_pentad(date_time=None):

    # compute dofyr and pentad for date_time

    date_time.dofyr = copy(date_time.day)
    # add up days in months prior to current month

    for i in arange(1,(date_time.month - 1)).reshape(-1):
        date_time.dofyr = copy(date_time.dofyr + days_in_month(date_time.year,i))

    date_time.pentad = copy(pentad_of_year(date_time.dofyr,date_time.year))
