using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using PhysicShapes;


public class two : MonoBehaviour
{
    public GameObject prefab;
    ThreeDBall ball;
    GameObject go;
    // Start is called before the first frame update
    void Start()
    {
        ball = new ThreeDBall(new Vector3(1, 0, .5f) * 3, Vector3.up * 8, .5f);
        go = Instantiate(prefab, ball.pos, Quaternion.identity);
    }

    // Update is called once per frame
    void FixedUpdate()
    {
        ball.simulate(Time.deltaTime * 3);
        go.GetComponent<Transform>().position = ball.pos;
    }
}
