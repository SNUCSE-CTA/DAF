/**
 * @file log.h
 * @author CTALAB
 * @brief log function macro
 * @version 0.1
 * @date 2021-03-29
 *
 * @copyright Copyright (c) 2021
 *
 */
#ifndef GLOBAL_LOG_H_
#define GLOBAL_LOG_H_

#define LOG_LEVEL 3

#ifdef PRINT_LOG
#define LOG(level, fmt, ...)                            \
  {                                                     \
    if (level >= LOG_LEVEL) printf(fmt, ##__VA_ARGS__); \
  }

#else

#define LOG(level, fmt, ...)

#endif  // PRINT_LOG
#endif  // GLOBAL_LOG_H_
