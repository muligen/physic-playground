using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class PhysicScene
{
    private static PhysicScene instance;
    private static readonly object lockObject = new object();
    private PhysicScene() {}
    public static PhysicScene Instance
    {
        get
        {
            // 使用双重检查锁定以确保线程安全
            if (instance == null)
            {
                lock (lockObject)
                {
                    if (instance == null)
                    {
                        instance = new PhysicScene();
                    }
                }
            }
            return instance;
        }
    }

    public float gravity = -9.8f;

    public float plane_lower_bound = 0f; 
    public float x_lower_bound = -10f; 
    public float x_upper_bound = 10f; 
    public float z_lower_bound = -10f; 
    public float z_upper_bound = 10f; 
}
