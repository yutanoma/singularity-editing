#pragma once

#include "./types.h"
#include "../store/geometry.h"
#include "../dispatcher/indices.h"

namespace rp {
namespace action {
// 選んだ2つのfaceのindexを交換する
IndexExchangeStateSetter use_index_exchange_setter(AngleSetter &angleSetter, IndicesSetter &indicesSetter);

IndexExchangeStateSetter use_add_singularity_pair_setter(AngleSetter &angleSetter, IndicesSetter &indicesSetter);

IndexExchangeStateSetter use_index_annihilation_setter(AngleSetter &angleSetter, IndicesSetter &indicesSetter);
}
}
