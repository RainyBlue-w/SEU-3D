require("/home/wuc/node_modules/core-js/actual/array/group-by.js")

const pick = (obj, arr) =>
    arr.reduce((iter, val) => (val in obj && iter.push(obj[val]), iter), [])

a = {
    'a': { text: 'A', label: 1 },
    'b': { text: 'B', label: 2 }
}

// console.log(pick(a, ['a', 'c']))

b = [
    { text: 'A', label: 1 },
    { text: 'B', label: 2 }
]

// console.log(b.map(obj => obj['label']))



const inventory = [
    { name: "芦笋", type: "蔬菜", quantity: 5 },
    { name: "香蕉", type: "水果", quantity: 0 },
    { name: "山羊", type: "肉", quantity: 23 },
    { name: "樱桃", type: "水果", quantity: 5 },
    { name: "鱼", type: "肉", quantity: 22 },
];

// var result = inventory.groupBy(({ type }) => type);

// console.log(result)

var result = inventory.groupBy( ({quantity}) => (quantity > 5 ? "ok" : "restock"))

console.log(result)


