# Usage:
# dump_data(metric, filename)
# where metric is one of:
#  * adjusted_moment_centroid
#  * area_per
# and more to come, and filename is the base of the filename of the shapefile
# (without the .shp/.dbf)
#
# TODO:
#  * clip to coastlines
#  * figure out if there's a sane thing to do with perimeter with
#    non-contractible shapes
#  * run on past and state data
#  * include population data, to compute population moments

import collections
import csv
import functools
import math
import shapefile
import os

def get_states(filename='data/state.txt'):
    """Returns a dict of FIPS codes to postal abbreviations."""
    with open(filename) as f:
        dialect = csv.Sniffer().sniff(f.read())
    with open(filename) as f:
        reader = csv.reader(f, dialect=dialect)
        next(reader) # ignore header
        return {state[0]: state[1] for state in reader}


def metric(function):
    """Decorator to make a function on pairs of points a function on shapes.

    function should take two points (each lists of length 2) and return some
    number, which will be added up over all consecutive pairs of points in the
    shape.

    Returns a function that takes a shape.  Intended to be used as a decorator.
    """
    @functools.wraps(function)
    def wrapped(shape):
        assert shape.shapeType == 5
        # Shapes consist of one or more parts, which begin at the positions
        # given in shape.parts.  The first part is the "outer" polygon, given
        # in clockwise order, and the latter parts are given in
        # counterclockwise order, and subtracted.  For most metrics, we really
        # just want to add up over all of them; perimeter is a bit wonky, so as
        # an approximation we simply ignore the removed bits.
        part_ends = shape.parts.tolist() + [len(shape.points)]
        pointss = [shape.points[part_ends[i]:part_ends[i+1]]
                   for i in range(len(shape.parts))]
        total = 0
        for points in pointss:
            # Note: points will actually start and end with the same point, but
            # we use it twice just in case.  If this screws up a metric, it's
            # probably a bad metric.
            prev = points[-1]
            for nxt in points:
                total += function(prev, nxt)
                prev = nxt
        return total
    return wrapped


@metric
def perimeter(prev, nxt):
    """Counts only the perimeter of the "outer" polygon, not any removed parts."""
    return math.sqrt((prev[0] - nxt[0]) ** 2 + (prev[1] - nxt[1]) ** 2)

def area_piece(prev, nxt):
    return nxt[0] * prev[1] - prev[0] * nxt[1]

@metric
def area(prev, nxt):
    """Area of a polygon with points in clockwise order."""
    return area_piece(prev, nxt) / 2

def area_per(points):
    """NGCH's Per_2: A/P^2, normalized to [0,1]."""
    return 4 * math.pi * area(points) / perimeter(points) ** 2

@metric
def center_x_helper(prev, nxt):
    return (prev[0] + nxt[0]) * area_piece(prev, nxt) / 6

@metric
def center_y_helper(prev, nxt):
    return (prev[1] + nxt[1]) * area_piece(prev, nxt) / 6

def centroid(points):
    a = area(points)
    return [center_x_helper(points)/a, center_y_helper(points)/a]

@metric
def moment_origin(prev, nxt):
    x_piece = (prev[0] * prev[0] + prev[0] * nxt[0] + nxt[0] * nxt[0])
    y_piece = (prev[1] * prev[1] + prev[1] * nxt[1] + nxt[1] * nxt[1])
    return (x_piece + y_piece) * area_piece(prev, nxt) / 12

def moment_centroid(points):
    centroid_x, centroid_y = centroid(points)
    return moment_origin(points) - area(points) * (
        centroid_x ** 2 + centroid_y ** 2)

def adjusted_moment_centroid(points):
    """NGCH's Dis_11, divided by sqrt(pi) to normalize it to be on [0,1]."""
    return area(points)/math.sqrt(2 * math.pi * moment_centroid(points))

def dump_data(metric, filebase='data/cb_2013_us_cd113_500k',
              states_filename='data/state.txt', out_filename=None):
    """Dumps a csv of state,district,metric.
    
    Notes:
      * States with at-large elections (i.e. a single rep) are denoted by
        district 00
      * Non-states with nonvoting reps (DC, PR, etc.) are listed and denoted by
        district 98 or 99
    """
    sf = shapefile.Reader(filebase)
    districts = sf.shapeRecords()
    states = get_states(states_filename)
    try:
        # data = [(dist.record[5], metric(dist.shape)) for dist in districts]
        data = [(states[dist.record[0]], dist.record[1], metric(dist.shape))
                for dist in districts]
    except:
        print dist.record
    data.sort()
    if out_filename is None:
        out_filename = filebase + '_' + metric.func_name + '.csv'
    with open(out_filename, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(data)

def block_map(state_fips, directory='/tmp/faces/'):
    all_files = os.listdir(directory)
    faces_files = [f for f in all_files
                   if f.startswith('tl_2014_%s' % state_fips)
                   and f.endswith('_faces.dbf')]
    if not faces_files:
        raise ValueError("No files found")
    b_map = collections.defaultdict(list)
    for f in faces_files:
        with open(directory + f, 'rb') as dbf:
            sf = shapefile.Reader(dbf=dbf)
            for rec in sf.records():
                block_code = ''.join(rec[1:4]) + rec[5]
                cd = rec[30]
                b_map[block_code].append(cd)
    cd_map = collections.defaultdict(list)
    for block_code in b_map:
        # Sometimes blocks get split; just take the district that most of the
        # pieces are in.  This isn't perfect but it's a small difference.
        district_options = collections.Counter(b_map[block_code])
        cd = district_options.most_common(1)[0][0]
        cd_map[cd].append(block_code)
    return cd_map

def block_reader(state_fips, directory='/tmp/faces/'):
    sf = shapefile.Reader(directory + 'tabblock2010_%s_pophu.shp' % state_fips)
    districts = enumerate(sf.records())
    indices_by_block_id = {record[4]: index for index, record in districts}
    def reader(block_id):
        if block_id in indices_by_block_id:
            return sf.shapeRecord(indices_by_block_id[block_id])
        else:
            print "block %s not found" % block_id
            # For some reason, a few blocks don't show up in the population
            # data.  Let's not worry about it.
            return None
    return reader

def population_centroid(block_ids, block_reader):
    centroid_x = 0
    centroid_y = 0
    total_pop = 0
    for block_id in block_ids:
        block = block_reader(block_id)
        if block:
            pop = block.record[7]
            block_centroid_x, block_centroid_y = centroid(block.shape)
            centroid_x += pop * block_centroid_x
            centroid_y += pop * block_centroid_y
            total_pop += pop
    return (centroid_x / total_pop, centroid_y / total_pop)

def population_moment(cd, block_ids, block_reader):
    moment = 0
    cd_centroid_x, cd_centroid_y = population_centroid(block_ids, block_reader)
    total_pop = 0
    for block_id in block_ids:
        block = block_reader(block_id)
        if block:
            pop = block.record[7]
            block_centroid_x, block_centroid_y = centroid(block.shape)
            centroid_dist = ((cd_centroid_x - block_centroid_x) ** 2 +
                             (cd_centroid_y - block_centroid_y) ** 2)
            moment += pop * centroid_dist
            total_pop += pop
    return area(cd.shape) * total_pop / (moment * 2 * math.pi)

def dump_population_moments(states=None,
                            filebase='data/tl_2014_us_cd114',
                            states_filename='data/state.txt',
                            block_data_directory='/tmp/faces/',
                            out_filename='data/population_moments.csv'):
    if states is None:
        # Filter out Puerto Rico etc. because there's no tabblock data for
        # them.
        states = [state for state in get_states().keys() if int(state) < 60]
    districts = shapefile.Reader(filebase).shapeRecords()
    # states -> dict of district number -> shapeRecord
    district_dict = collections.defaultdict(dict)
    for district in districts:
        district_dict[district.record[0]][district.record[1]] = district
    for state in states:
        data = []
        print "processing %s" % get_states()[state]
        bm = block_map(state)
        br = block_reader(state)
        for num in district_dict[state]:
            try:
                data.append((get_states()[state], num,
                             population_moment(district_dict[state][num],
                                               bm[num], br)))
            except Exception as e:
                print "%s error in %s-%s" % (e, get_states()[state], num)
        data.sort()
        with open(out_filename, 'a') as f:
            writer = csv.writer(f)
            writer.writerows(data)


