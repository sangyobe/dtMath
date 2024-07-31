#include "testPrint.h"

void PrintTitle(const char* str)
{
    char asterisk[81];
    char title[81];
    memset(asterisk, '*', sizeof(asterisk));
    memset(title, ' ', sizeof(title));

    asterisk[0] = '/';
    asterisk[78] = '/';
    asterisk[79] = '\n';
    asterisk[80] = 0;

    title[0] = '/';
    title[1] = '*';
    title[77] = '*';
    title[78] = '/';
    title[79] = '\n';
    title[80] = 0;

    memcpy(&title[3], str, strlen(str));

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
    Printf("\033[92m");
    Printf("%s", asterisk);
    Printf("%s", title);
    Printf("%s", asterisk);
    Printf("\033[0m");
#else
    Printf("%s", asterisk);
    Printf("%s", title);
    Printf("%s", asterisk);
#endif
}

void PrintHeading(const char* str)
{
    char title[81];
    memset(title, '*', sizeof(title));

    title[0] = '/';
    title[1] = '*';
    title[2] = ' ';
    title[78] = '/';
    title[79] = '\n';
    title[80] = 0;
    memcpy(&title[3], str, strlen(str));

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
    Printf("\033[36m");
    Printf("%s", title);
    Printf("\033[0m");
#else
    Printf("%s", title);
#endif
}