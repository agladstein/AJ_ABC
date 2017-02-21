from bisect import bisect_right

def find2(a, x):
    """This function receives the array with available sites (sites that passed the frequency cut-off)"""
    if x <= a[0]: #x is smaller than first value in a
        return 0
    elif x >= a[len(a)-1]: #x is larger than last value in a
        return len(a)-1
    #Find leftmost item greater than or equal to x
    else:
        i = bisect_right(a, x)
        j = i-1
        d1=abs(a[i]-x)
        d2=abs(a[j]-x)
        if d2<d1:
            return j
        elif d1<d2:
            return i
        elif d1==d2:
            return i

def add_snps(avail_sites, nb_avail_sites, pos_asc, nbss_asc, nb_array_snps):
    first_index=pos_asc[0]
    last_index=pos_asc[-1]
    if(nb_avail_sites>nbss_asc): ###this should happen all the time
        if (last_index<nb_avail_sites-1):

            try:
                avail_sites[last_index+1]
            except:
                print "well, it WASN'T defined after all!"
                try:
                    avail_sites[first_index-1]
                except:
                    print "well, it WASN'T defined after all!"
                else:
                    for i in xrange(len(pos_asc)):
                        pos_asc[i]=pos_asc[i]-1
            else:
                print "sure, it was defined"
                pos_asc.append(last_index+1)

        elif (last_index==nb_avail_sites-1):

            if(first_index-1)>=0:
                pos_asc.insert(0,first_index-1)
    return pos_asc
