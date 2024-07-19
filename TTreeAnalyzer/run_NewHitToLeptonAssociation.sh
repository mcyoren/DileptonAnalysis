#!/bin/sh

root -l -b << EOF
    .L MyEvent_OLD.C+
    .L ../AnaTrain/Run14AuAuLeptonComby/MyEvent.C+
    .x NewHitAssociation.C
EOF