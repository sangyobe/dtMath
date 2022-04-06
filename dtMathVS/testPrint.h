#pragma once

#include "./dtMath/dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdio.h>
#include <string.h>
#define Printf      printf
#define Println     printf("\n")
#elif defined(ARDUINO)
#include <Arduino.h>
#define Printf      Serial.printf
#define Println     Serial.println()
#endif

void PrintTitle(const char* str);
void PrintHeading(const char* str);