#include "testGnuPlot.h"
#include "./dtMath/dtGnuPlot.hpp"

void Test_GnuPlot()
{
    CdtGnuPlot<> plot1;
    CdtGnuPlot<> plot2;

    const int dataLen = 500;
    float time = 0;
    float dataT[dataLen];
    float data1[dataLen];
    CdtVector<dataLen> data2;

    for (int i = 0; i < dataLen; i++)
    {
        dataT[i] = time;
        data1[i] = sin(time);
        data2(i) = cos(time);
        time += 0.01f;
    }

    plot1.Cmd("plot [-pi/2:pi] cos(x),-(sin(x) > sin(x+1) ? sin(x) : sin(x+1))");

    plot2.SetXlabel("time", 13, "Times");
    plot2.SetYlabel("amp", 13, "Times");
    plot2.SetTitle("Sin and Cos Curve", 15, "Times");
    plot2.LinePoint(dataT, data1, dataLen, "cos", "0x0000FF", 1, plot2.PT_CIRCLE, 1, -5);
    plot2.Dash(dataT, data2.GetElementsAddr(), dataLen, "sin", "-", "red", 1);
    plot2.Draw();
}