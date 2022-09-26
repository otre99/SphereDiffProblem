
def sec_to_strtime(secs):
    if secs < 60.0:
        return "{:.1f} secs".format(secs)
    elif secs < 60*60.0:
        return "{:.2f} mins".format(secs/60.0)
    elif secs < 60.0*60.0*24.0:
        return "{:.2f} hrs".format(secs/60.0/60.0) 
    else:
        return "{:.2f} days".format(secs/60.0/60.0/24.0)                 

def sec_to_days(secs):
    return secs/(24*60*60)