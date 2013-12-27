import sys, operator

from random import random

import bisect

verbose = False
output = sys.stdout

class ParseError(Exception):
    def __init__(self, value="Parsing Error"):
        self.parameter = value
    def __str__(self):
        return self.parameter

def prefix_filename(filename, prefix):
    """Add a prefix to a filename. 

    i.e. filename = ./test_data/anns.bed, prefix = complement
    returns ./test_data/complement_anns.bed
    """
    ## TODO I should be using the sys path handling machinery here
    ## XXX This is probably broken on windows
    import re
    filename_portions = re.split("/", filename)
    filename_portions[-1] = prefix + filename_portions[-1]
    return '/'.join( filename_portions )

##################### Python2.3 Compatability Code #############################
try: set
except NameError:
    from sets import Set as set
##################### END Python2.3 Compatability Code ##########################

##################### Standard Statistical Methods #############################
### Functions we need if numpy isn't installed
###
################################################################################

try: from numpy.random import multivariate_normal
except ImportError:
    def multivariate_normal(meanMat, covMat, n):
        """ this constructs 2D multivariate dist, given a mean and cov list 
        """

        from math import sqrt
        from random import gauss as normal

        # make sure that the mean and cov entries make sense
        assert len(meanMat) == 2
        assert len(covMat) == 2
        assert len(covMat[0]) == 2
        assert len(covMat[1]) == 2

        # make sure that cov matrix is symmetric
        assert covMat[0][1] == covMat[1][0]

        # first, calculate the cholesky decomposition of the cov matrix
        # I'll use an explicit Cholesky-Crout algorithm since I know it's 2x2
        L11 = sqrt( covMat[0][0] )
        L12 = covMat[0][1]/L11
        # I add +0.000001 to covMat[1][1] - L12**2 because rounding errors sometimes
        # give me very small negative values, and then sqrt raises a domain error
        # I could use a try-catch, but I like to avoid the overhead since this
        # is a potentially hot spot in the algorithm - and it doesnt really matter
        L22 = sqrt( covMat[1][1] - L12**2 + 0.000001)

        mvn = [ [], [] ]
        # then I'll build a 2xn vector of MVN random samples
        for loop in xrange(n):
            z1 = normal(0,1)
            z2 = normal(0,1)

            mvn[0].append( meanMat[0] + L11*z1 )
            mvn[1].append( meanMat[1] + L12*z1 + L22*z2 )

        return zip(mvn[0], mvn[1])

try: from scipy.stats.stats import mean, cov, std
except ImportError:
    from math import sqrt
    def mean(data): return float(sum(data))/len(data)

    def var(data):
        ev = mean(data)
        return sum( [ (entry - ev)**2 for entry in data ] )/( len(data) - 1 )

    def std(data): return sqrt( var(data) )

    def cov(data):
        """ Takes an iterable of 'vectors' ( by which I mean iterables of the same length )
        """
        # make sure they are all the same length
        for entry in data: assert len(entry) == len(data[0])

        # its expecting to get the entries in the opposite order that we want them
        data = zip(*data)

        # pairwise covariance
        def _cov(X,Y):
            assert len(X) == len(Y)
            mx = mean(X)
            my = mean(Y)
            return sum( [ float((x-mx)*(y-my)) for x,y in zip(X,Y) ] )/( len(X) - 1 )

        # first, make a zeros covariance matrix
        covMat = [ [0]*len(data) for loop in xrange(len(data)) ]
        # fill in the covaraince matrix entries
        for oloop in xrange(len(data)):
            for iloop in xrange(oloop, len(data)):
                tmpVar = _cov(data[oloop], data[iloop])
                covMat[iloop][oloop] = tmpVar
                covMat[oloop][iloop] = tmpVar
        return covMat

try:
    from scipy.stats.distributions import norm as norm_dist
    sn_cdf = norm_dist(loc=0, scale=1).cdf
except ImportError:
    def sn_cdf(x):
        # brute force integration of the sn pdf
        # probably accumualtes errors and has all sorts of other
        # undesirable properties
        # however, it is still accurate to 1e-8 for the worst case
        #
        # this could be much faster/more accurate if I weighted
        # the points better

        # if we're more than 20 sd over the mean, return 0.0
        if x > 20: return 0.0

        # if we're more than 20 sd less than the mean, return 0
        if x < -20: return 1.0

        def norm_pdf(x, u, s):
            from math import pi, sqrt, exp

            constTerm = 1/(s*sqrt(2*pi))
            expTerm = exp(-0.5*(s**-2)*((x-u)**2))

            return constTerm*expTerm

        # if we are past 20, the added value to the cdf
        # is unsubstantial: so just integrate to 20
        a = max(-20, float(x))
        b = 20.0

        estimate = 0
        num = 10000
        for loop in xrange(num):
            ln = a + loop*((b-a)/num)
            hn = a + (loop+1)*((b-a)/num)
            local_est = (hn-ln)*norm_pdf((hn+ln)/2.0, 0, 1)
            estimate += local_est

        # this integrates above x - so reverse it
        estimate = 1-estimate

        ## some code to test my approximation
        #from scipy.stats.distributions import norm as norm_dist
        #sci_sn_cdf = norm_dist(loc=0, scale=1).cdf
        #print estimate, sci_sn_cdf(x), sci_sn_cdf(x) - estimate

        return estimate

##################### End Standard Statistical Methods #########################


########################### cumulative_data Type ###############################

class cumulative_data(object):
    def __init__(self, data, length=None, typ=None):
        """This initializes a cummulative array. 

        """ 
        self.split_points = None
        ## Note that the start_feature_cumsu should be len(feature_regions)+1
        ## this is asserted later on the code
        # store the feature regions
        self.feature_regions = []
        # store the cumsum at the bp *previous* to the region start
        self.start_feature_cumsum = []

        ## counters to store the current index and cummulative sum
        currSum = 0
        index = 0

        ## set the length to be the region length
        self.length = data.length

        ## check and see if data implements the iter_feature_regions method
        ## if it does, then we will only store the feature boundaries. This
        ## is a memory optimization for binary feature sets
        if hasattr(data, 'iter_feature_regions'):
            for feature_interval in data.iter_feature_regions():    
                self.feature_regions.append(feature_interval)
                self.start_feature_cumsum.append(currSum)
                currSum += feature_interval.size
            ## add a final cumulative sum for the end points
            self.start_feature_cumsum.append(currSum)

        ## otherwise, store every single feature. This is *slow* for big data
        ## sets. So slow in fact, that we will raise a warning if we expect this
        ## to take over (arbitrarily) 100 megs of memory TODO
        else:
            for entry in data.iter_features():
                self.feature_regions.append(interval(index, index))
                self.start_feature_cumsum.append(currSum)
                currSum += entry
                index += 1
            ## add a final cumulative sum for the end points
            self.start_feature_cumsum.append(currSum)

    def __len__(self):
        return self.length

    def iter_split_points(self, minValue, maxValue):
        iterator = iter(self.feature_regions)
        interval = iterator.next()
        while interval.end < minValue:  
            interval = iterator.next()
        while interval.start < maxValue:  
            yield interval.start
            if interval.size > 1 and interval.end < maxValue:
                yield interval.end
            interval = iterator.next()                
        return

    def __getitem__(self, index):
        # make sure that we're only trying to index by basepair
        assert isinstance(index, (int, long))
        l_index = bisect.bisect_left(self.feature_regions, index)

        ## if the l_index is greater than any regions, we want the full cumsum
        if l_index == len(self.feature_regions):
            return self.start_feature_cumsum[l_index]
        ## else, find the start of the next interval
        else: 
            next_interval_start = self.feature_regions[l_index].start
        ## if we are in a no region area, return the cumsum at the next point
        if next_interval_start > index:
            return self.start_feature_cumsum[l_index]
        ## if we are in a region area, return the cumsum at the next point
        ## plus 1 for every bp inside the region we are in. Note that we add one
        ## to this because the intervals are inclusive
        else:
            return self.start_feature_cumsum[l_index] + (index - next_interval_start + 1)

    def plot(self, split_points=[]):
        import pylab

        feature_starts = [ interval.start for interval in self.feature_regions ]
        feature_starts.append(self.length)
        pylab.plot(feature_starts, self.start_feature_cumsum)
        for entry in split_points:
            pylab.axvline(entry)
        pylab.show()

######################### END cumulative_data Type #############################

########################### pointwise data fns #################################

class pointwise_data(list):
    """Store pointwise data.

    This is meant to store the data from a pointwise file. It is just a list of
    tuples 

    TO-DO: override append and set-item methods to make sure the data stays 
    sorted. ( Is this really necessary/helpful? )

    Methods:
        add - add a new index a data point.
        shift - make all of the indexes consecutive
        find_common_indexes - find the indexes that this pw object shares w/
            another pw object
        update_from_indexes - copy this pw object including only the indexes
            in the passed argument
        to_file - write a pw file from the data stored in this file
        iter_features - iterate through all of the score data


    """
    def add(self, index, score):
        #assert ( len(self) == 0 or ( self[-1][0] <= index ) )
        self.append((index, score))


    def shift(self):
        """Shift the indexes to start at zero and be consecutive.

        """
        rv = pointwise_data()

        for index in xrange(len(self)):
            rv.add(index, self[index][1])

        return rv

    def find_common_indexes(self, other):
        """Returns a sorted list of indexes that exist in both self and other

        """
    
        #assert isinstance(other, pointwise_data)
        # BAD python2.3 compat change
        selfInd = set([ index for index, score in self ])
        # BAD python2.3 compat change
        otherInd = set([ index for index, score in other ])

        rv = list( selfInd.intersection(otherInd) )
        rv.sort()

        return rv

    def update_from_indexes(self, indexes):
        """Returns a copy of pointwise_data containing only the values that 
           have indexes in the passed list.

        """

        assert isinstance(indexes, list)

        # sort the indexes for the merge
        indexes.sort()

        rv = pointwise_data()

        loop = 0
        for index in indexes:
            while( self[loop][0] < index ):
                loop += 1

            if self[loop][0] == index:
                rv.add( index, self[loop][1] )

        return rv

    def to_file(self, filename):
        """Write the data in this object to a pointwise file.

        """
        of = open(filename, 'w')

        for index, score in self:
            of.write("%i\t%f\n" % (index, score))
        
        of.close()

    def iter_features(self):
        """Iterate through the scores.

        """
        for index, score in self:
            yield score

def parse_pointwise_file(fp):
    """Parse a pointwise file into a pointwise_data object.

    """

    import re

    data = pointwise_data()

    for line in fp:
        index, score = re.split('\s+', line.strip())[0:2]
        data.add(int(index), float(score))
        
    return data

def generate_pointwise_data(length, indexExistsProb=0.5, valueMean=2, fp=sys.stdout):
    """Generate random pointwise data.

    """

    import random
    for index in xrange(int(length/indexExistsProb)):
        if random.random() > indexExistsProb:
            fp.write("%i\t%f\n" % ( index, 2*valueMean*random.random() ))

########################### END pointwise data fns #############################

class interval(tuple):
    """Store a closed interval of ints. 
    
    The interval is closed at both ends. ie, interval(1,3) is {1,2,3}
    """
    # this prevents each object from creating a __dict__ ( just a small optimization )
    __slots__ = []

    # creates the class from the object
    def __new__(self, start, end):
        if end < 0 or start < 0 or  not isinstance(start, (int, long)) != int or not isinstance(end, (int, long)):
            raise ValueError, 'All values must be non-negative (integers|longs), not (%i, %i)' % ( start, end )
        if end < start:
            raise ValueError, 'The end value (%d) must be greater than the start value (%d) ' % ( end, start )
        self = tuple.__new__(interval, (start, end))
        return self

    start = property(lambda self: self[0])
    end = property(lambda self: self[1])
    size = property(lambda self: self[1] - self[0]+1 ) 

    # calculate the length of the overlap
    def overlap(self, otherInterval):
        assert isinstance(otherInterval, interval)
        return max(min(self[1], otherInterval[1]) - max(self[0], otherInterval[0]) + 1, 0)

    def does_overlap(self, otherInterval):
        assert isinstance(otherInterval, interval)
        return self.overlap(otherInterval) > 0

    ############### Comparators ########################################################################################
    def __lt__(self, other):
        if isinstance(other, (int, long, float)):
            return self.end < other
        elif type(other) == interval:
            return self.end < other.start
        else: raise TypeError, "Can't compare an object of type %s to a slice object" % type(other)

    def __le__(self, other):
        if isinstance(other, (int, long, float)):
            return self.start <= other
        elif type(other) == interval:
            return self.end <= other.end
        else: raise TypeError, "Can't compare an object of type %s to a slice object" % type(other)

    def __eq__(self, other):
        # only allow the comparison of intervals or 'numbers'
        if not isinstance(other, (int, long, float, interval)): return False
        return self <= other and self >= other

    def __ne__(self, other):
        return not self == other

    def __ge__(self, other):
        if isinstance(other, (int, long, float)):
            return self.end >= other
        elif type(other) == interval:
            return self.start >= other.start
        else: raise TypeError, "Can't compare an object of type %s to a slice object" % type(other)

    def __gt__(self, other):
        if isinstance(other, (int, long, float)):
            return self.start > other 
        elif type(other) == interval:
            return self.start > other.end
        else: raise TypeError, "Can't compare an object of type %s to a slice object" % type(other)

class region(list):
    def __init__(self, iterator=(), name=None):
        list.__init__(self, iterator)
        self.name = name
    
    # override the append method so that it keeps the elements sorted
    def append(self, newItem):
        # if the item belongs at the end, insert it
        if len(self) > 0 and newItem > self[-1]:
            list.append(self, newItem)
        # if the item doesn't belong at the end make sure it's unique and
        # insert it into the appropriate location
        else:
            # first, make sure that bisect_right == bisect_left
            r_loc = bisect.bisect_right(self, newItem)            
            l_loc = bisect.bisect_left(self, newItem)            
            # if the new element intersects a previous element
            if r_loc != l_loc:
                new_interval = interval( min(self[l_loc].start, newItem.start), max(self[r_loc-1].end, newItem.end) )
                if verbose and False:
                    print "WARNING! WARNING! WARNING!\n \tFeature interval %s intersects %i other FI(s) " % ( newItem, r_loc-l_loc )
                    print "\tThe intersecting intervals will be merged."
                    print "\t%s is being changed to %s" % ( self[l_loc], new_interval )
                for index in xrange(r_loc-1, l_loc, -1):
                    if verbose and False:
                        print "\tThe interval %s is being removed." % str(self[index])
                    del self[index]
                if verbose and False: print "WARNING! WARNING! WARNING!\n"
                self[l_loc] = new_interval
                #raise ValueError, 'The new element %s intersects a previous element' % str(newItem)
            else: list.insert(self, r_loc, newItem)

    # same as append: actually add is probably the most appropriate
    def insert(self, newItem): self.append(newItem)
    def add(self, newItem): self.append(newItem)

    def numRegions(self):
        return len(self)

    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop, None))

    def __getitem__(self, key):
        ## if the key is an integer, treat it like a normal list
        if type(key) == int:
            return list.__getitem__(self, key)        
        
        ## otherwise, the key better be a slice
        assert type(key) == slice
        return type(self)( iter(list.__getitem__(self, key)), self.name )

class genomic_coverage_region( region ):
    def _length(self):
        try: return self._cached_length
        except AttributeError:
            # BAD python 2.3 compat change
            self._cached_length = sum( [ feature.size for feature in self ] )
        return self._cached_length
    length = property( _length )  

    def append(self, newItem):
        # clear the cached length 
        if hasattr(self, "_cached_length"):
            del self._cached_length
        # call region's append method 
        region.append(self, newItem)

    def intersection(self, inter):
        """Return the intersection of inter and self, *in genomic coordinates*
        
        inter needs to be an interval.
        """

        # first, find the lower and upper bounds of intersection
        l_loc = bisect.bisect_left(self, inter.start)
        r_loc = bisect.bisect_right(self, inter.end)        
        # if there is no overlap, return None
        if l_loc == r_loc: return None

        # make sure that self has *some* regions
        assert len(self) > 0
        # make sure the interval is before the first basepair
        assert self[0].start <= inter.end
        # makre sure the interval is after the last basepair
        assert self[-1].end >= inter.start

        # otherwise, calculate the overlap
        # BAD python 2.3 compat change
        bp_overlap = sum( [item.overlap(inter) for item in self] )
        # next, note the interval start is the count of bp's <= self[l_loc] 
        if l_loc > 0:
            new_start = (self[0:l_loc]).length
        else:
            new_start = 0
        # account for the potential overlap
        if inter.start > self[l_loc].start:
            new_start += ( inter.start - self[l_loc].start )
        
        # make sure that the new start is positive - ( explicitly catch ian's bug )
        assert new_start >= 0

        return interval( new_start, new_start+bp_overlap-1 )

class feature_region( region ):
    """Store the features for a specific region.

    """

    def __init__(self, length, iterator, name=None):
        #for interval in iterator:
        #    self.add(interval)
        list.__init__(self, iterator)
        self.length = length
        self.name = name

    def randomSubSample(self, percLen):
        startIndex = int(round(random()*self.length))
        endIndex = int(round(startIndex + percLen*self.length))
        return self[startIndex:endIndex]

    # calculates the overlap between this and another set of intervals track
    def overlap(self, otherAnnTrack):
        # special case for empty intervals        
        if len(self) == 0 or len(otherAnnTrack) == 0: return 0

        totalOverlap = 0

        currentMatches = feature_region(max(self.length, otherAnnTrack.length), ())
        
        thisIter = iter(self)
        nextMatch = thisIter.next()
        
        for iv in otherAnnTrack:
            # first, get rid of any current matches that dont match
            # because of the ordering, these should start not matching at the begining
            for item in currentMatches:
                if item.end >= iv.start: break
                else: del currentMatches[0] 

            # next, add any new items to currentMatches
            while nextMatch != None and nextMatch.start <= iv.end:
                currentMatches.append(nextMatch)
                try: nextMatch = thisIter.next()
                except StopIteration: nextMatch = None

            # finally, calculate the overlap of every item in currentMatches
            for item in currentMatches:
                totalOverlap += item.overlap(iv)

        return totalOverlap

    def featuresLength(self):
        # BAD python 2.3 compat change
        return sum( [ feature.size for feature in self ] )

    # calculates the region overlap between this and another set of intervals track
    def regionOverlap(self, otherRegion):
        """Calculate the number of regions that overlap self with the given region

        We take self to be the covered region. ie, for 
        self    ===============================
        other   ---   ----   -----------    -         
        the region overlap is 1, whereas for  
        self    ===========  ==================
        other   ---   ----   -----------    -                 
        it is 2.

        The algorithm is nearly identical to they overlap fn - it's a merge join
        """

        # special case for empty intervals        
        if len(self) == 0: return 0

        totalOverlap = 0
        self_index = 0
        other_index = 0
         
        for iv in self:
            # remove the other intervals that can never overlap (given the sort)
            while other_index < len(otherRegion) and \
                  otherRegion[other_index].end < iv.start:
                other_index += 1
            # if other index is longer than the length, there can be no matches
            if other_index >= len(otherRegion):
                break

            # check to see if there is an overlap
            if iv.does_overlap( otherRegion[other_index] ):
                totalOverlap += 1
                # if so, move onto the next region in self
                continue
            # if not, we know that other is *not* << self and there is not
            # any overlap, so we continue tot he next self anyways.

        #if self.brute_regionOverlap( otherRegion ) != totalOverlap:
        #    print "WARNING!!! NOT EQUAL!!!!", self.brute_regionOverlap( otherRegion ), totalOverlap

        return totalOverlap

    def bruteOverlap(self, otherRegion):
        # mostly for testing my fancier overlap algorithm
        # this compares every possible combination of intervals
        totalOverlap = 0
        for other_iv in otherRegion:
            for self_iv in self:
                totalOverlap += other_iv.overlap(self_iv)
        return totalOverlap

    def brute_regionOverlap(self, otherRegion):
        # mostly for testing my fancier overlap algorithm
        # this compares every possible combination of regions
        totalOverlap = 0
        for self_iv in self:
            for other_iv in otherRegion:
                if other_iv.does_overlap(self_iv):
                    totalOverlap += 1
                    break
        return totalOverlap

    def iter_features(self):
        """Iterates through every feature in the interval.
        
        For instance, if the region was of length 10 and had a single
        feature interval [2,8] this would return 0 0 1 1 1 1 1 1 1 0.
        """

        # deal with the start points
        for loop in xrange(self[0][0]):
            yield 0

        # deal with the internal features
        for loop in xrange(len(self)-1):
            for iloop in xrange(self[loop][1] - self[loop][0]):
                yield 1
            for iloop in xrange(self[loop+1][0] - self[loop][1]):
                yield 0

        # deal with the last feature
        for loop in xrange(self[-1][1] - self[-1][0]):
            yield 1

        # deal with the non-features until the end
        for loop in xrange(self.length - self[-1][1]):
            yield 0

        return

    def iter_feature_regions(self):
        """Iterates through every feature region in the interval.
        
        For instance, if the region was of length 10 and had a single
        feature interval [2,8] this would return interval(2, 8).
        """

        for feature_interval in self:
            yield feature_interval

        return

    def split(self, split_points):
        """Splits this region into a regions object.

        """
        rv = regions()
        
        # deal with the trivial case of no splits
        if len(split_points) == 0:
            rv[self.name] = self
            return rv

        name = ""
        if self.name != None:
            name = self.name + "_"

        rv[name + 'split_%i' % 0] = feature_region(split_points[0], self[0:split_points[0]], name + 'split_%i' % 0)

        for index in xrange(1, len(split_points)):
            rv[name + 'split_%i' % index] = feature_region(
                split_points[index]-split_points[index-1], 
                self[(split_points[index-1]):(split_points[index]+1)], 
                name + 'split_%i' % index)

        rv[name + 'split_%i' % len(split_points)] = feature_region(
                self.length - split_points[len(split_points)-1],
                self[(split_points[-1]):(self.length+1)], 
                name + 'split_%i' % len(split_points))

        assert self.length == sum( item.length for item in rv.values() )

        return rv

    def __getslice__(self, start, stop):
        return self.__getitem__(slice(start, stop, None))

    def __getitem__(self, key):
        ## if the key is an integer, treat it like a normal list
        if type(key) == int:
            return list.__getitem__(self, key)        
        
        ## otherwise, the key better be a slice
        assert type(key) == slice
        start = key.start
        stop = key.stop
        ## if the start or stop are less than 1, then interpret it as a percentage
        if type(start) == float and start < 1:
            start = int(start*self.length)
        if type(stop) == float and stop <= 1:
            stop = int(stop*self.length)

        if start == None: start = 0
        if stop == None: stop = self[-1].end

        l_loc = bisect.bisect_left(self, start)
        r_loc = bisect.bisect_right(self, stop)

        ## take the list slice
        ## note that this is explicitly modifying the end points
        # BAD python 2.3 compat change
        gen = [ interval(max(instart-start,0), min(instop-start, stop)) for instart, instop in list.__getitem__(self, slice(l_loc, r_loc)) ]
        retValue = feature_region(stop-start, gen) 

        return retValue

from random import uniform
class regions(dict):
    """A region container object. 
    
    This is just a dictionary of region's whereby each key is the 
    region's track name. Also contains a method to store an overlap stat.
    """

    #!!! BUG - FIX THE LENGTH CACHING MECHANISM
    # It currently assumes that I never use __setitem__

    def _totalLength(self):
        """Calculate the total length and cache it for later use.

        """
        try: return self._cached_totalLength
        except AttributeError:
            self._cached_totalLength = float(sum([ value.length for value in self.values() ]))
            return self._cached_totalLength

    totalLength = property( _totalLength )

    def regionFraction(self):
        """Returns a dict of the relative region lengths.

        For instance, if there are two regions of names R1 and R2 and they are both
        the same length then this will return { 'R1': 0.5, 'R2': 0.5 }

        """
        # BAD python2.3 compat change
        return dict([ ( key, float(value.length)/self.totalLength ) for key, value in self.iteritems() ])

    def extend(self, other):
        for key in other.keys():
            # make sure that we wont overwrite an existing region
            if self.has_key(key):
                raise ValueError, 'Can not extend %r with %r: they both contain a region named %s' % (self, other, key)
            self[key] = other[key]

        # delete the cached length
        try: 
            del self._cached_totalLength
        except AttributeError:
            pass

    def writeBedFile(self, bed_f, lengths_f):
        keys = self.keys()
        keys.sort()

        for key in keys:
            for start, end in self[key]:
                bed_f.write("%s\t%i\t%i\n" % (key, start, end))

            lengths_f.write("%s\t0\t%i\n" % (key, self[key].length) )


def parse_bed_file(bed_f, lengths_file):
    """Parse a bed and lengths file into a regions object.

    """

    import re

    # first, build the lengths file region
    lengths = regions()
    for line_num, line in enumerate(lengths_file):
        values = re.split('\W+', line.strip())
        if len(values) < 3:
            raise ParseError("Error parsing line #%i \"%s\" in the file \"%s\"" % (line_num+1, line.strip(), lengths_file.name))  
        chName, start, end = values[0:3]
            
        feature = interval(int(start), int(end))
        # if the lengths object doesnt have a region of this name, create it
        if not lengths.has_key(chName):
            lengths[chName] = genomic_coverage_region(name=chName)
        lengths[chName].append(feature)

    # next, initialize the chromosomes data structure
    chromosomes = regions()
    for chName in lengths.keys():
        length = lengths[chName].length
        chromosomes[chName] = feature_region(length, (), chName)

    for line in bed_f:
        # parse out the region name, and interval and parse it
        chName, start, end = re.split('\W+', line.strip())[0:3]
        feature = interval(int(start), int(end))
        
        # make sure that a genomic region with name chName exists in chromosomes
        # if it doesn't, raise a warning and continue - we implcitly assume that 
        # region not being in lengths means we dont care about it
        if not lengths.has_key(chName):
            print "WARNING: The named region '%s' is in the bed file but NOT the domain file." % chName
            print "WARNING: This could indicate that the domain file is incorrect."
            continue

        # see what parts - if any - of the feature intersects the genomic region
        # take only the intersecting pieces and shift the coordinates to be
        # measured in the genomic region space
        shifted_feature = lengths[chName].intersection(feature) 
        # if there is no intersection, there is nothing else to do
        if shifted_feature == None: continue
        else: chromosomes[chName].append(shifted_feature)

    lengths_file.seek(0)
    return chromosomes
