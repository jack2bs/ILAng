#include <stdio.h>
#include <string.h>

int main(void)
{
    int runningSum = 0;
    char strBuf[128];
    char discardStr[128];
    int numCycles;
    int discardInt;
    double avgFUDynPower;
    double avgMemDynPower;
    double avgFULeakPower;
    double avgMemLeakPower;
    double fuArea;
    double memArea;
    double discardDouble;

    FILE * outFile = fopen("ppa_models.csv", "w");

    fprintf(outFile, "func,cycles,avg fu dynamic power,fu leakage power,avg mem dynamic power,mem leakage power\n");

    while (1)
    {
        runningSum = 0;
        memset(strBuf,0,strlen(strBuf));
        
        runningSum += scanf("===============================\n");
        runningSum += scanf("        Aladdin Results        \n");
        runningSum += scanf("===============================\n");
        runningSum += scanf("Running : %127s\n", discardStr);
        runningSum += scanf("Top level function: %127s\n", strBuf);
        runningSum += scanf("Cycle : %d cycles\n", &numCycles);
        runningSum += scanf("Upsampled Cycle : %d cycles\n", &discardInt);
        runningSum += scanf("Avg Power: %lf mW\n", &discardDouble);
        runningSum += scanf("Idle FU Cycles: %d cycles\n", &discardInt);
        runningSum += scanf("Avg FU Power: %lf mW\n", &discardDouble);
        runningSum += scanf("Avg FU Dynamic Power: %lf mW\n", &avgFUDynPower);
        runningSum += scanf("Avg FU leakage Power: %lf mW\n", &avgFULeakPower);
        runningSum += scanf("Avg MEM Power: %lf mW\n", &discardDouble);
        runningSum += scanf("Avg MEM Dynamic Power: %lf mW\n", &avgMemDynPower);
        runningSum += scanf("Avg MEM Leakage Power: %lf mW\n", &avgMemLeakPower);
        runningSum += scanf("Total Area: %lf uM^2\n", &discardDouble);
        runningSum += scanf("FU Area: %lf uM^2\n", &fuArea);
        runningSum += scanf("MEM Area: %lf uM^2\n", &memArea);
        runningSum += scanf("%*[^=]");
        runningSum += scanf("===============================\n");
        runningSum += scanf("        Aladdin Results        \n");
        runningSum += scanf("===============================\n");

        if (runningSum <= 0)
        {
            break;
        }
        

        fprintf(outFile, "%s,%d,%lf,%lf,%lf,%lf\n", strBuf, numCycles, avgFUDynPower, avgFULeakPower, avgMemDynPower, avgMemLeakPower);     


    }
    fclose(outFile);

}