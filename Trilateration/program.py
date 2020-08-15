import math
import numpy
import argparse
import csv


# python main.py --input input.csv --output outputCsvFile (kullanim)


def trilateration(dist_line):
    yaricapDunya = 6371
    DistA = float(dist_line[1]) / 1000
    DistB = float(dist_line[2]) / 1000
    DistC = float(dist_line[3]) / 1000
    LatA = float(dist_line[5])
    LonA = float(dist_line[4])
    LatB = float(dist_line[7])
    LonB = float(dist_line[6])
    LatC = float(dist_line[9])
    LonC = float(dist_line[8])

    # Otantik küre (authalic sphere :)

    #  Enlem ve Boylam verilerini ECEF Koordinat Sistemine uygun şekilde çevirme

    #   1. Enlem ve boylamın radyana çevrilmesi
    #   2. Radyan verilerinin ECEF koordinat sistemine çevrilmesi

    xA = yaricapDunya * (math.cos(math.radians(LatA))
                         * math.cos(math.radians(LonA)))
    yA = yaricapDunya * (math.cos(math.radians(LatA))
                         * math.sin(math.radians(LonA)))
    zA = yaricapDunya * (math.sin(math.radians(LatA)))

    xB = yaricapDunya * (math.cos(math.radians(LatB))
                         * math.cos(math.radians(LonB)))
    yB = yaricapDunya * (math.cos(math.radians(LatB))
                         * math.sin(math.radians(LonB)))
    zB = yaricapDunya * (math.sin(math.radians(LatB)))

    xC = yaricapDunya * (math.cos(math.radians(LatC))
                         * math.cos(math.radians(LonC)))
    yC = yaricapDunya * (math.cos(math.radians(LatC))
                         * math.sin(math.radians(LonC)))
    zC = yaricapDunya * (math.sin(math.radians(LatC)))

    P1 = numpy.array([xA, yA, zA])
    P2 = numpy.array([xB, yB, zB])
    P3 = numpy.array([xC, yC, zC])

    # başlangıç noktasındaki daire için dönüşüm
    # x eksenin daire2 elde edilmesi için gerekli olan dönüşüm

    ex = (P2 - P1) / (numpy.linalg.norm(P2 - P1))
    i = numpy.dot(ex, P3 - P1)
    ey = (P3 - P1 - i * ex) / (numpy.linalg.norm(P3 - P1 - i * ex))
    ez = numpy.cross(ex, ey)
    d = numpy.linalg.norm(P2 - P1)
    j = numpy.dot(ey, P3 - P1)

    # değerlerin kullanılması
    x = (math.pow(DistA, 2) - math.pow(DistB, 2) + math.pow(d, 2)) / (2 * d)
    y = ((math.pow(DistA, 2) - math.pow(DistC, 2)
          + math.pow(i, 2) + math.pow(j, 2)) / (2 * j)) - ((i / j) * x)

    try:
        z = math.sqrt(math.pow(DistA, 2) - math.pow(x, 2) - math.pow(y, 2))
    except:
        z = float('nan')

    # triPt nesnesi üçleme noktasının ECEF x,y,z' ye sahip bir dizisidir.
    triPt = P1 + x * ex + y * ey + z * ez

    # Enlem boylam derecelerinin ECEF'e çevrilmesi ve ardından dereceye dönüştürülmesi

    enlem = math.degrees(math.asin(triPt[2] / yaricapDunya))
    boylam = math.degrees(math.atan2(triPt[1], triPt[0]))

    return [dist_line[0], boylam, enlem]


#  Komut satırı üzerinden csv çıktısı elde etmek için gerekli olan tanımlamalar
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i', '--input', help='Input file',
                               required=True)
    requiredNamed.add_argument('-o', '--output', help='Output TSV file',
                               required=True)
    args = parser.parse_args()

    with open(args.output, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['Data', 'Longitude', 'Latitude'])
        with open(args.input, 'r') as f:
            for line in f:
                out_line = trilateration(line.strip().split(','))
                writer.writerow(out_line)
