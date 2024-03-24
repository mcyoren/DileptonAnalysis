#ifndef __DECLARATIONS_H__
#define __DECLARATIONS_H__

#define SQR(x) ((x) * (x))
#define central(x) ((x) < 0 ? 0 : (x) > 4 ? 4  : (x))
#define REG_HIST(dim, h, ...) \
    h = new TH##dim##D(#h, #h, __VA_ARGS__); \
    se->registerHisto(#h, h)

#define REG_HISTOS(dim, h, n, ...) \
    for (int i = 0; i < n; ++i) { \
        char buffer[128] = {0}; \
        sprintf(buffer, #h "_%d", i); \
        h[i] = new TH##dim##D(buffer, buffer, __VA_ARGS__); \
        se->registerHisto(buffer, h[i]); \
    }

#define REG_HISTOS_2N(dim, h, n, k, ...) \
    for (int i = 0; i < n; ++i) { \
        for (int j = 0; j < k; j++) \
        { \
            char buffer[128] = {0}; \
            sprintf(buffer, #h "_%d_%d", i, j); \
            h[i][j] = new TH##dim##D(buffer, buffer, __VA_ARGS__); \
            se->registerHisto(buffer, h[i][j]); \
        } \
    }

#endif
