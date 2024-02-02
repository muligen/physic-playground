using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class one : MonoBehaviour
{
    public Image Image;
    public const float gravity = -9.8f;
    public Vector3 speed = new Vector3(10, 0, 0);
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void FixedUpdate()
    {
        var pos = Image.GetComponent<Transform>().position;
        speed += new Vector3(0, Time.deltaTime * gravity, 0);
        if(pos.x < 0 || pos.x > 1920f)
        {
            speed.x = -speed.x;
        }
        if (pos.y < 0 || pos.y > 1080f)
        {
            speed.y = -speed.y;
        }
        Image.GetComponent<Transform>().position += speed;
    }
}
