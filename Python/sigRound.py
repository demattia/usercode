# #!/usr/bin/python

# Evaluates the digits to be used with the round function in order to have
# the first two significant figures


def sigRound(error):
    tempErr = error
    digits = 0
    if( int(tempErr) == 0 ):
        while (int(tempErr) == 0):
            tempErr = tempErr*10
            digits += 1
    else:
        digits = 1
        while (int(tempErr) != 0):
            tempErr = tempErr/10
            digits -= 1

    # Check for rounding
    if( int(round(error,digits)*10**(digits-1)) != 0 ):
        return digits
    return digits+1

# Round the number to the two significant figures of the error
##a = 0.1235131223183
##aErr = 0.0012392882843

##digits = sigRound(aErr)
##print "digits =", digits

##a = round(a,digits)
##aErr = round(aErr,digits)

##print a, "+-", aErr
