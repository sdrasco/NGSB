#LIBS = ../lib/libNGSB.a /usr/lib/librecipes_c.a
LIBS = ../lib/libNGSB.a ../lib/librecipes_c.a
CC = gcc
#CFLAGS = -g -pg  # debug & profile
#CFLAGS = -O2 -pg  # optimize the speed & profile
CFLAGS = -O2  # optimize the speed
#CFLAGS = -O3  # optimize the hell out of the speed
#INCLUDES = -I../include/ -I/usr/local/recipes_c-ansi/include
#INCLUDES = -I../include/  -I/usr/local/src/recipes_c-ansi/include
INCLUDES = -I../include/

all: BurstDistance.o BurstStatistic.o GetAlpha.o Power.o PairApproxLogLikelihood.o PairApproxLogLikelihood2.o MaxLikelihood.o NonGaussianDistance.o GaussianDistance.o CrossCorr.o GaussianStatistic.o Gaussian.o NGSBData.o PairNGSBData.o ExactLogLikelihood.o PairExactLogLikelihood.o SimplexSort.o SimplexSort2.o BoundaryCheck.o BoundaryCheck2.o SimplexMaxLikelihood.o SimplexMaxLikelihood2.o SimplexMaxLikelihood3.o GeneralLogLikelihood.o Approx1.o Approx2.o SeedRand.o SimplexTest SimplexTest2 GaussianTest NGSBDataTest DetectMeasure MultiDetect FileMultiDetect LogLogMultiDetect PairLogLogMultiDetect ApproxLogLogMultiDetect GeneralMultiDetect StandardPerformance ComparePerformance ComparePerformance2 ComparePerformance3 DetectionCurves ExactApprox timing ErrorScale BurstTest FastDetectionCurves BurstDistribution

BurstDistance.o: BurstDistance.c
	$(CC) $(CFLAGS) $(INCLUDES) -c BurstDistance.c -o BurstDistance.o;
	ar -vr ../lib/libNGSB.a BurstDistance.o

BurstStatistic.o: BurstStatistic.c
	$(CC) $(CFLAGS) $(INCLUDES) -c BurstStatistic.c -o BurstStatistic.o;
	ar -vr ../lib/libNGSB.a BurstStatistic.o

GetAlpha.o: GetAlpha.c
	$(CC) $(CFLAGS) $(INCLUDES) -c GetAlpha.c -o GetAlpha.o;
	ar -vr ../lib/libNGSB.a GetAlpha.o

Power.o: Power.c
	$(CC) $(CFLAGS) $(INCLUDES) -c Power.c -o Power.o;
	ar -vr ../lib/libNGSB.a Power.o

PairApproxLogLikelihood.o: PairApproxLogLikelihood.c
	$(CC) $(CFLAGS) $(INCLUDES) -c PairApproxLogLikelihood.c -o PairApproxLogLikelihood.o;
	ar -vr ../lib/libNGSB.a PairApproxLogLikelihood.o

PairApproxLogLikelihood2.o: PairApproxLogLikelihood2.c
	$(CC) $(CFLAGS) $(INCLUDES) -c PairApproxLogLikelihood2.c -o PairApproxLogLikelihood2.o;
	ar -vr ../lib/libNGSB.a PairApproxLogLikelihood2.o

MaxLikelihood.o: MaxLikelihood.c
	$(CC) $(CFLAGS) $(INCLUDES) -c MaxLikelihood.c -o MaxLikelihood.o;
	ar -vr ../lib/libNGSB.a MaxLikelihood.o

CrossCorr.o: CrossCorr.c
	$(CC) $(CFLAGS) $(INCLUDES) -c CrossCorr.c -o CrossCorr.o;
	ar -vr ../lib/libNGSB.a CrossCorr.o

NonGaussianDistance.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c NonGaussianDistance.c -o NonGaussianDistance.o;
	ar -vr ../lib/libNGSB.a NonGaussianDistance.o

GaussianDistance.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c GaussianDistance.c -o GaussianDistance.o;
	ar -vr ../lib/libNGSB.a GaussianDistance.o

GaussianStatistic.o: GaussianStatistic.c
	$(CC) $(CFLAGS) $(INCLUDES) -c GaussianStatistic.c -o GaussianStatistic.o;
	ar -vr ../lib/libNGSB.a GaussianStatistic.o

Gaussian.o: Gaussian.c
	$(CC) $(CFLAGS) $(INCLUDES) -c Gaussian.c -o Gaussian.o;
	ar -vr ../lib/libNGSB.a Gaussian.o

NGSBData.o: NGSBData.c
	$(CC) $(CFLAGS) $(INCLUDES) -c NGSBData.c -o NGSBData.o;
	ar -vr ../lib/libNGSB.a NGSBData.o

PairNGSBData.o: PairNGSBData.c
	$(CC) $(CFLAGS) $(INCLUDES) -c PairNGSBData.c -o PairNGSBData.o;
	ar -vr ../lib/libNGSB.a PairNGSBData.o

ExactLogLikelihood.o: ExactLogLikelihood.c
	$(CC) $(CFLAGS) $(INCLUDES) -c ExactLogLikelihood.c -o ExactLogLikelihood.o;
	ar -vr ../lib/libNGSB.a ExactLogLikelihood.o

PairExactLogLikelihood.o: PairExactLogLikelihood.c
	$(CC) $(CFLAGS) $(INCLUDES) -c PairExactLogLikelihood.c -o PairExactLogLikelihood.o;
	ar -vr ../lib/libNGSB.a PairExactLogLikelihood.o

SimplexMaxLikelihood.o: SimplexMaxLikelihood.c
	$(CC) $(CFLAGS) $(INCLUDES) -c SimplexMaxLikelihood.c -o SimplexMaxLikelihood.o;
	ar -vr ../lib/libNGSB.a SimplexMaxLikelihood.o;

SimplexMaxLikelihood2.o: SimplexMaxLikelihood2.c
	$(CC) $(CFLAGS) $(INCLUDES) -c SimplexMaxLikelihood2.c -o SimplexMaxLikelihood2.o;
	ar -vr ../lib/libNGSB.a SimplexMaxLikelihood2.o;

SimplexMaxLikelihood3.o: SimplexMaxLikelihood3.c
	$(CC) $(CFLAGS) $(INCLUDES) -c SimplexMaxLikelihood3.c -o SimplexMaxLikelihood3.o;
	ar -vr ../lib/libNGSB.a SimplexMaxLikelihood3.o;

SimplexSort.o: SimplexSort.c
	$(CC) $(CFLAGS) $(INCLUDES) -c SimplexSort.c -o SimplexSort.o;
	ar -vr ../lib/libNGSB.a SimplexSort.o

SimplexSort2.o: SimplexSort2.c
	$(CC) $(CFLAGS) $(INCLUDES) -c SimplexSort2.c -o SimplexSort2.o;
	ar -vr ../lib/libNGSB.a SimplexSort2.o

BoundaryCheck.o: BoundaryCheck.c
	$(CC) $(CFLAGS) $(INCLUDES) -c BoundaryCheck.c -o BoundaryCheck.o;
	ar -vr ../lib/libNGSB.a BoundaryCheck.o

BoundaryCheck2.o: BoundaryCheck2.c
	$(CC) $(CFLAGS) $(INCLUDES) -c BoundaryCheck2.c -o BoundaryCheck2.o;
	ar -vr ../lib/libNGSB.a BoundaryCheck2.o

GeneralLogLikelihood.o: GeneralLogLikelihood.c
	$(CC) $(CFLAGS) $(INCLUDES) -c GeneralLogLikelihood.c -o GeneralLogLikelihood.o;
	ar -vr ../lib/libNGSB.a GeneralLogLikelihood.o

Approx1.o: Approx1.c
	$(CC) $(CFLAGS) $(INCLUDES) -c Approx1.c -o Approx1.o;
	ar -vr ../lib/libNGSB.a Approx1.o

Approx2.o: Approx2.c
	$(CC) $(CFLAGS) $(INCLUDES) -c Approx2.c -o Approx2.o;
	ar -vr ../lib/libNGSB.a Approx2.o

SeedRand.o: SeedRand.c
	$(CC) $(CFLAGS) $(INCLUDES) -c SeedRand.c -o SeedRand.o;
	ar -vr ../lib/libNGSB.a SeedRand.o

Krand.o: Krand.c
	$(CC) $(CFLAGS) $(INCLUDES) -c Krand.c -o Krand.o;
	ar -vr ../lib/libNGSB.a Krand.o

SimplexTest: SimplexTest.c
	$(CC) $(CFLAGS) SimplexTest.c -o SimplexTest $(INCLUDES) $(LIBS) -lm

SimplexTest2: SimplexTest2.c
	$(CC) $(CFLAGS) SimplexTest2.c -o SimplexTest2 $(INCLUDES) $(LIBS) -lm

GaussianTest: GaussianTest.c
	$(CC) $(CFLAGS) GaussianTest.c -o GaussianTest $(INCLUDES) $(LIBS) -lm

NGSBDataTest: NGSBDataTest.c
	$(CC) $(CFLAGS) NGSBDataTest.c -o NGSBDataTest $(INCLUDES) $(LIBS) -lm

DetectMeasure: DetectMeasure.c
	$(CC) $(CFLAGS) DetectMeasure.c -o DetectMeasure $(INCLUDES) $(LIBS) -lm

MultiDetect: MultiDetect.c
	$(CC) $(CFLAGS) MultiDetect.c -o MultiDetect $(INCLUDES) $(LIBS) -lm

FileMultiDetect: FileMultiDetect.c
	$(CC) $(CFLAGS) FileMultiDetect.c -o FileMultiDetect $(INCLUDES) $(LIBS) -lm

LogLogMultiDetect: LogLogMultiDetect.c
	$(CC) $(CFLAGS) LogLogMultiDetect.c -o LogLogMultiDetect $(INCLUDES) $(LIBS) -lm

PairLogLogMultiDetect: PairLogLogMultiDetect.c
	$(CC) $(CFLAGS) PairLogLogMultiDetect.c -o PairLogLogMultiDetect $(INCLUDES) $(LIBS) -lm

ApproxLogLogMultiDetect: ApproxLogLogMultiDetect.c
	$(CC) $(CFLAGS) ApproxLogLogMultiDetect.c -o ApproxLogLogMultiDetect $(INCLUDES) $(LIBS) -lm

GeneralMultiDetect: GeneralMultiDetect.c
	$(CC) $(CFLAGS) GeneralMultiDetect.c -o GeneralMultiDetect $(INCLUDES) $(LIBS) -lm

StandardPerformance: StandardPerformance.c
	$(CC) $(CFLAGS) StandardPerformance.c -o StandardPerformance $(INCLUDES) $(LIBS) -lm

ComparePerformance: ComparePerformance.c
	$(CC) $(CFLAGS) ComparePerformance.c -o ComparePerformance $(INCLUDES) $(LIBS) -lm

ComparePerformance2: ComparePerformance2.c
	$(CC) $(CFLAGS) ComparePerformance2.c -o ComparePerformance2 $(INCLUDES) $(LIBS) -lm

ComparePerformance3: ComparePerformance3.c
	$(CC) $(CFLAGS) ComparePerformance3.c -o ComparePerformance3 $(INCLUDES) $(LIBS) -lm

DetectionCurves: DetectionCurves.c
	$(CC) $(CFLAGS) DetectionCurves.c -o DetectionCurves $(INCLUDES) $(LIBS) -lm

ExactApprox: ExactApprox.c
	$(CC) $(CFLAGS) ExactApprox.c -o ExactApprox $(INCLUDES) $(LIBS) -lm

timing: timing.c
	$(CC) $(CFLAGS) timing.c -o timing $(INCLUDES) $(LIBS) -lm

ErrorScale: ErrorScale.c
	$(CC) $(CFLAGS) ErrorScale.c -o ErrorScale $(INCLUDES) $(LIBS) -lm

BurstTest: BurstTest.c
	$(CC) $(CFLAGS) BurstTest.c -o BurstTest $(INCLUDES) $(LIBS) -lm

FastDetectionCurves: FastDetectionCurves.c
	$(CC) $(CFLAGS) FastDetectionCurves.c -o FastDetectionCurves $(INCLUDES) $(LIBS) -lm

BurstDistribution: BurstDistribution.c
	$(CC) $(CFLAGS) BurstDistribution.c -o BurstDistribution $(INCLUDES) $(LIBS) -lm

