// debug.h
#pragma once

// Method 1: enable via compile option -DDEBUG
#ifdef DEBUG
#define DEBUG_LOG(x)                                                           \
  do {                                                                         \
    std::cout << "[DEBUG] " << x << std::endl;                                 \
  } while (0)
#else
#define DEBUG_LOG(x)                                                           \
  do {                                                                         \
  } while (0) // No-op; the compiler will optimize it away
#endif

// Method 2: custom switch (recommended)
#ifndef ENABLE_DEBUG_LOG
#define ENABLE_DEBUG_LOG 0 // Disabled by default
#endif

#if ENABLE_DEBUG_LOG
#define DEBUG_LOG(x)                                                           \
  do {                                                                         \
    std::cerr << "[DEBUG][" << __FILE__ << ":" << __LINE__ << "] " << x        \
              << std::endl;                                                    \
  } while (0)
#else
#define DEBUG_LOG(x)                                                           \
  do {                                                                         \
  } while (0)
#endif
